"""
CLI commands
------------

    python comp.py [model1] [model2]

Run the full `model1` and `model2`

 ---

    python comp.py model1 mech_sweep

Repeat N `model1` runs incrementally inserting mechanisms, where N is the
total numbers of mechanisms.

 ---

    python comp.py model1 skip:mech1

Run a mechanism sweep but skip all runs that don't insert `mech1`, consider
the following series:

mechA
mechA mechB
mechA mechB mechC

 ---

    python comp.py model1 mech_sweep skip:mechB

would skip ahead past the first run introducing `mechB` to the `mechA mechB mech C` run

 ---

    python comp.py model1 mech_sweep blocklist:mechA;mechB;mechC

Run a mechanism sweep, filtering out `mechA`, `mechB` and `mechC`
"""

import dbbs_models, nrnsub
import arbor, sys, itertools
import plotly.graph_objs as go
from arborize import compose_types, flatten_composite
import arborize.core

def monkey_apply_section_ions(self, section, ions):
    prop = {"e": "e{}", "int": "{}i", "ext": "{}e"}
    for ion_name, ion_props in ions.items():
        for prop_name, value in ion_props.items():
            try:
                prop_attr = prop[prop_name].format(ion_name)
            except KeyError as e:
                raise IonAttributeError(f"Unknown ion attribute '{prop_name}'.") from None
            try:
                setattr(section.__neuron__(), prop_attr, value)
            except:
                pass

arborize.core.NeuronModel._apply_section_ions = monkey_apply_section_ions

class single_recipe(arbor.recipe):
    def __init__(self, cell, probes, Vm, K):
        # (4.1) The base C++ class constructor must be called first, to ensure that
        # all memory in the C++ class is initialized correctly.
        arbor.recipe.__init__(self)
        self.the_cell = cell
        self.the_probes = probes
        self.the_props = arbor.neuron_cable_properties()
        self.the_props.set_property(Vm=Vm, tempK=K, rL=35.4, cm=0.01)
        self.the_props.set_ion(ion='na', int_con=10,   ext_con=140, rev_pot=50)
        self.the_props.set_ion(ion='k',  int_con=54.4, ext_con=2.5, rev_pot=-77)
        self.the_props.set_ion(ion='ca', int_con=0.00005, ext_con=2, rev_pot=132.5)
        self.the_props.set_ion(ion='h', valence=1, int_con=1.0, ext_con=1.0, rev_pot=-34)

        self.the_cat = arbor.default_catalogue()
        self.the_cat.extend(arbor.dbbs_catalogue(), "")
        self.the_props.register(self.the_cat)

    def num_cells(self):
        # (4.2) Override the num_cells method
        return 1

    def num_sources(self, gid):
        # (4.3) Override the num_sources method
        return 0

    def cell_kind(self, gid):
        # (4.4) Override the cell_kind method
        return arbor.cell_kind.cable

    def cell_description(self, gid):
        # (4.5) Override the cell_description method
        return self.the_cell

    def probes(self, gid):
        # (4.6) Override the probes method
        return self.the_probes

    def global_properties(self, kind):
        # (4.7) Override the global_properties method
        return self.the_props


def whitelist_model(model, whitelist):
    class model(model): pass
    model.section_types = model.section_types.copy()
    for k, v in model.section_types.items():
        model.section_types[k] = v = v.copy()
        try:
            mech_iter = v["mechanisms"].items()
        except:
            v._self["mechanisms"] = {mech_key: mech for mech_key, mech in mech_iter if mech_key in whitelist}
        else:
            v["mechanisms"] = {mech_key: mech for mech_key, mech in mech_iter if mech_key in whitelist}

    return model

def run_nrn(setup, model, duration, v_init, temperature):
    from patch import p
    neuron_cell = model()
    neuron_cell.record_soma()
    neuron_time = p.time
    p.dt = 0.025
    p.cvode.active(0)
    p.celsius = temperature - 273.15
    p.finitialize(v_init)
    p.continuerun(duration)
    return neuron_time, neuron_cell.Vm


def run_arb(setup, model, duration, v_init, temperature):
    cell = model.cable_cell(Vm=v_init)
    probe = arbor.cable_probe_membrane_voltage('(proximal (tag 1))')
    recipe = single_recipe(cell, [probe], v_init, temperature)
    context = arbor.context()
    domains = arbor.partition_load_balance(recipe, context)
    sim = arbor.simulation(recipe, domains, context)
    sim.record(arbor.spike_recording.all)
    handle = sim.sample((0, 0), arbor.regular_schedule(0.1))
    sim.run(tfinal=duration)

    spikes = sim.spikes()
    data, meta = sim.samples(handle)[0]
    return data[:, 0], data[:, 1]


def run_model(model, duration=1000, v_init=-65, setup=None, temperature=305.15):
    if setup is not None:
        if "whitelist" in setup:
            whitelist = setup["whitelist"]
            model = whitelist_model(model, whitelist)
    return run_nrn(setup, model, duration, v_init, temperature), run_arb(setup, model, duration, v_init, temperature)

if __name__ == "__main__":
    models = {
        v: {"model": getattr(dbbs_models, v)}
        for v in sys.argv
        if v.endswith("Cell")
    } or {
        name: {"model": model}
        for name, model in vars(dbbs_models).items()
        if name.endswith("Cell")
    }
    blocklister = [v for v in sys.argv if v.startswith("blocklist:")]
    if blocklister:
        blocklist = set("".join(blocklister[0].split(":")[1:]).split(";"))
    else:
        blocklist = set()
    for name, setup in models.items():
        model = setup["model"]
        if "mech_sweep" in sys.argv:
            mechs = sorted(list(set(itertools.chain(*(flatten_composite(model, v).get("mechanisms", dict()).keys() for v in model.section_types.values())))), key=str)
            skippers = [k.split(":")[1] for k in sys.argv if "skip:" in k]
            skip = list(str(m) for m in mechs).index(skippers[0]) + 2 if skippers else 0
            setup["setups"] = [{"whitelist": [m for m in mechs[:i] if str(m) not in blocklist]} for i in range(skip, len(mechs))]
        setups = setup.get("setups", (None,))
        for setup in setups:
            if setup is not None:
                sname = name + " " + ", ".join(map(str, setup["whitelist"]))
            else:
                sname = name
            print("Running", sname)
            (nt, nv), (at, av) = run_model(model, setup=setup)
            print("Setup", sname, "finished")
            go.Figure(
                [
                    go.Scatter(
                        x=nt,
                        y=nv,
                        name=f"NEURON: {sname}"
                    ),
                    go.Scatter(
                        x=at,
                        y=av,
                        name=f"Arbor: {sname}"
                    ),
                ],
                layout_title_text=sname,
            ).show()
            if "one" in sys.argv:
                exit()
