from probeinterface import ProbeGroup, write_prb

def make_probe_from_recording(rec,probe_path):
    probe = rec.get_probe() 
    pg = ProbeGroup()
    pg.add_probe(probe)
    write_prb(probe_path, pg)
    print(probe)