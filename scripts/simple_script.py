import sys, os

sys.path.append("../")
try:
    from scrubber.core import ScrubberCore
except Exception as e:
    print("The script must be executed from the script/ directory")
    raise e
# sys.exit(1)


if __name__ == "__main__":
    import sys
    from pprint import pprint as pp

    config = ScrubberCore.get_defaults()
    # print("CONFIG!")
    with open("default_config.py", "w") as fp:
        pp(config, fp)
    # pp(config)

    config["input"]["values"]["fname"] = sys.argv[1]
    # config["errors"]["values"]["input_err_fname"] = "errors_proc.input"
    # config["errors"]["values"]["process_err_fname"] = "errors_proc.sdf"
    # config["errors"]["values"]["process_err_ftype"] = "sdf"
    # the protein
    # config["input"]["values"]["name_property"] = "NON_EXISTENT_PROPERTY"
    # config["input"]["values"]["start_count"] = 70000
    # config["input"]["values"]["start_count"] = 160
    config["input"]["values"]["end_count"] = 50

    config["general"]["values"]["max_proc"] = 8
    config["general"]["values"]["nice"] = 40

    config["isomers"]["active"] = True
    config["isomers"]["values"]["verbose"] = False
    config["isomers"]["values"]["proto_enum"] = True
    config["isomers"]["values"]["proto_keep_all"] = False
    config["isomers"]["values"]["proto_pH"] = 7.4  # single value
    config["isomers"]["values"]["ph_datafile"] = "..//scrubber/data/test_model.txt"
    # config["isomers"]["values"]["proto_pH"] = [6.4, 8.4]  # pH range

    # config["isomers"]["values"]['stereoisomer_enum'] = 'undefined' #  'all', False
    config["isomers"]["values"]["stereo_enum"] = False
    # 'undefined', 'all', 'False'
    config["isomers"]["values"]["tauto_enum"] = True

    # # test neutralizing already charged molecules
    # config["isomers"]["values"]['proto_enum'] = True
    # config["isomers"]["values"]['proto_neutralize'] = True

    # reject molecules that did not converge
    config["geometry"]["values"]["strict"] = True

    # common output format settings
    config["output"]["values"]["fname"] = "dest/testing.sdf"
    config["output"]["values"]["fname"] = "tautomers_test/tautomers.sdf"
    config["output"]["values"]["fname"] = "tautomers.sdf"
    config["output"]["values"]["ftype"] = "sdf"
    # single SDF file mode
    if False:
        config["output"]["values"]["mode"] = "single"
    else:
        # multi SDF mode
        config["output"]["values"]["mode"] = "split"
        config["output"]["values"]["disable_name_sanitize"] = True
        config["output"]["values"]["naming"] = "name"
        config["output"]["values"]["fname"] = "out/"
        config["output"]["values"]["ftype"] = "sdf"
        # config["output"]["values"]["naming_field"] = "MISSING"
        config["output"]["values"]["max_lig_per_dir"] = 100

    # sys.exit()
    sc = ScrubberCore(options=config)
    sc.process_file()
    # graph =  sc.isomer.reaction_log
    # traj = graph.trajectory
    # from pprint import pprint as pp
    # pp(traj)
    # graph.print_history()


### TNOTES
# in single mode, the filename must be specified
# un split mode no
