
- Change path inside the codes:
    1. in cross_link.py:
        OBJ_DIR - directory that contains the xl_objects pickle files
        OBJ_DIR_PREFIX - prefix of OBJ_DIR
        PROCESSED_OBJ_DIR - directory that contains xl_objects pivkle with distance info
        PROCESSED_OBJ_DIR_PREFIX - prefix of PROCESSED_OBJ_DIR
    2. in pdb_files_manger.py:
        PDB_FILES_DIR - directory that will contain pdb files
        INTER_AND_INTRA_LYS_XL_FILES_PATH - path to directory that needed in data generation process
        INTER_INTRA_LYS_XL_NEIGHBORS_FILES_PATH - path to directory that needed in data generation process
        XL_NEIGHBORS_FEATURE_DICT_INTER_INTRA_LYS - path to feature dictionary. needed in data generation process.
        XL_NEIGHBORS_EXE - path to the xl_neighbors executable (C++ program that extracts the features
    3. in general_utils.py:
        OBJ_DIR - directory that contains the pickle files
        PROJ_DIR - the project dir
    4. in alpha_fold_files.py:
        AF_PDB_DIR - directory for alpha fold pdb files
        AF_PRED_ERROR_DIR - directory for alpha fold pae files

- Generate the dataset on your own from csv file:
    1. In order to create objects of crosslinks from csv:
        a. If your csv file is from XlinkNet github / in the same format:
            If the csv file contains distances (it does in github) run:
                python3 cross_link.py --run_mode import --input_file_path <path_to_csv> --final_save_name <name_of_res_file> --download_files <True/False>
            Else:
                python3 cross_link.py --run_mode import --is_processed False --out_name <path_of_saving_processed> --input_file_path <path_to_csv> --final_save_name <name_of_res_file> --download_files <True/False>
        b. If your input file is in other format, please check the code in general_xl_parser.py and run:
                python3 cross_link.py --run_mode parse --parse_cfg <path_to_config_file> --is_processed False --out_name <path_of_saving_processed> --input_file_path <path_to_csv> --final_save_name <name_of_res_file> --download_files <True/False>

    2. Extract features from cross link objects (takes few hours, you can download features / dataset from github)
        python3 cross_link.py --run_mode extract --final_save_name <name_of_saved_object_list>

    3. Generate pytorch dataset for NN - run:
        python3 train.py --create_dataset True --xl_objects_file <name_of_saved_object_list> --dataset_name <name_for_the_new_dataset>

- Now that the dataset is set, we can train the model by running the train script along with the config file that we want:
    python3 train.py --cfg <path_to_config_file>