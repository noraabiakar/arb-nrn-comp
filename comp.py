import dbbs_models as models
import arbor, sys
from patch import p
import plotly.graph_objs as go
import time

class single_recipe(arbor.recipe):
    def __init__(self, cell, probes):
        # (4.1) The base C++ class constructor must be called first, to ensure that
        # all memory in the C++ class is initialized correctly.
        arbor.recipe.__init__(self)
        self.the_cell = cell
        self.the_probes = probes
        self.the_props = arbor.neuron_cable_properties()
        self.the_props.set_property(Vm=-65, tempK=300, rL=35.4, cm=0.01)
        self.the_props.set_ion(ion='na', int_con=10,   ext_con=140, rev_pot=50)
        self.the_props.set_ion(ion='k',  int_con=54.4, ext_con=2.5, rev_pot=-77)
        self.the_props.set_ion(ion='ca', int_con=0.0001, ext_con=2, rev_pot=132.5)
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


def run_nrn(setup, model, duration, v_init):
    if setup is not None:
        pass
    neuron_cell = model()
    neuron_cell.record_soma()
    neuron_time = p.time
    p.dt = 0.025
    p.cvode.active(0)
    p.celsius = 32
    p.finitialize(-40)
    p.continuerun(1000)
    return neuron_time, neuron_cell.Vm


def run_arb(setup, model, duration, v_init):
    if setup is not None:
        pass
    cell = model.cable_cell()
    probe = arbor.cable_probe_membrane_voltage('(proximal (tag 1))')
    recipe = single_recipe(cell, [probe])
    context = arbor.context()
    domains = arbor.partition_load_balance(recipe, context)
    sim = arbor.simulation(recipe, domains, context)
    sim.record(arbor.spike_recording.all)
    handle = sim.sample((0, 0), arbor.regular_schedule(0.025))
    sim.run(tfinal=duration)

    spikes = sim.spikes()
    data, meta = sim.samples(handle)[0]
    return data[:, 0], data[:, 1]


def run_single_model(model, duration=1000, v_init=-40):
    print("run_nrn")
    start = time.time()
    nrn_result = run_nrn(None, model, duration, v_init)
    end = time.time()
    print("time: ", end - start, "s")

    print("run_arb")
    start = time.time()
    arb_result = run_arb(None, model, duration, v_init)
    end = time.time()
    print("time: ", end - start, "s")

    return nrn_result, arb_result
    #return run_nrn(None, model, duration, v_init), run_arb(None, model, duration, v_init)

def run_model_setups(setups, model, duration=1000, v_init=-40):
    result = {}
    for k, setup in setups.items(): 
        print("run_nrn")
        nrn_result = run_nrn(setup, model, duration, v_init)
        print("arb_nrn")
        arb_result = run_arb(setup, model, duration, v_init)
        print("done")
        result[k] = (nrn_result, arb_result) 
    return results
    #return {k: (run_nrn(setup, model, duration, v_init), run_arb(setup, model, duration, v_init)) for k, setup in setups.items()}

if __name__ == "__main__":
    setups = {
        v: {"model": getattr(models, v)}
        for v in sys.argv
        if v.endswith("Cell")
    } or {
        name: {"model": model}
        for name, model in vars(models).items()
        if name.endswith("Cell")
    }
    for name, setup in setups.items():
        model = setup["model"]
        print("Running", name)
        (nt, nv), (at, av) = run_single_model(model)
        print("Setup", name, "finished")
        go.Figure(
            [
                go.Scatter(
                    x=nt,
                    y=nv,
                    name=f"{name} - NEURON"
                ),
                go.Scatter(
                    x=at,
                    y=av,
                    name=f"{name} - Arbor"
                ),
            ]
        ).show()
