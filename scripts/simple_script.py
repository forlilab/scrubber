import sys, os
sys.path.append("../")
try:
    from scrubber.core.scrubbercore import ScrubberCore
except Exception as e:
    print("The script must be executed from the script/ directory")
    raise e
# sys.exit(1)


if __name__ == "__main__":
    import sys
    from pprint import pprint as pp

    config = ScrubberCore.get_defaults()
    print("CONFIG!")
    pp(config)

    config["input"]["values"]["fname"] = sys.argv[1]
    # the protein
    # config["input"]["values"]["name_property"] = "NON_EXISTENT_PROPERTY"
    # config["input"]["values"]["start_count"] = 70000
    # config["input"]["values"]["start_count"] = 160
    # config["input"]["values"]["end_count"] = 166


    config["general"]["max_proc"] = 8
    config["general"]["nice"] = 40

    config["isomers"]["active"] = True
    config["isomers"]["values"]['verbose'] = False
    config["isomers"]["values"]['ph_datafile'] = "..//scrubber/data/test_model.txt"
    config["isomers"]["values"]['protomer_keep_all'] = False
    config["isomers"]["values"]['protomer_enum'] = False
    config["isomers"]["values"]['protomer_pH'] = 7.4        # single value
    config["isomers"]["values"]['protomer_pH'] = [6.4, 8.4] # pH range

    # config["isomers"]["values"]['stereoisomer_enum'] = 'undefined' #  'all', False
    config["isomers"]["values"]['stereoisomer_enum'] = False # 'undefined', 'all', 'False'
    config["isomers"]["values"]['tautomer_enum'] = True

    # # test neutralizing already charged molecules
    # config["isomers"]["values"]['protomer_enum'] = True
    # config["isomers"]["values"]['protomer_neutralize'] = True

    # reject molecules that did not converge
    config["geometry"]["values"]['strict'] = True

    # common output format settings
    config["output"]["values"]["out_fname"] = "dest/testing.sdf"
    config["output"]["values"]["out_fname"] = "tautomers_test/tautomers.sdf"
    config["output"]["values"]["out_fname"] = "tautomers.sdf"
    config["output"]["values"]["out_format"] = "sdf"
    # single SDF file mode
    if False:
        config["output"]["values"]["mode"] = "single"
    else:
        # multi SDF mode
        config["output"]["values"]["mode"] = "split"
        config["output"]["values"]["sanitize_name"] = True
        config["output"]["values"]["naming"] = "name"
        config["output"]["values"]["out_fname"] = "out/"
        config["output"]["values"]["out_format"] = "sdf"
        # config["output"]["values"]["naming_field"] = "MISSING"
        config["output"]["values"]["max_lig_per_dir"] = 100

    # sys.exit()
    sc = ScrubberCore(options=config)
    sc.process()
    # graph =  sc.isomer.reaction_log
    # traj = graph.trajectory
    # from pprint import pprint as pp
    # pp(traj)
    # graph.print_history()



### TNOTES
# in single mode, the filename must be specified
# un split mode no

