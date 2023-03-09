import numpy as np
import general_utils
import yaml
import cross_link
import csv
import data_table_parser


def load_config(cfg_file):
    with open(cfg_file, "r") as f:
        cfg = yaml.safe_load(f)
        return cfg


class XlParser:
    def __init__(self, cfg):
        self.deliminator = cfg['DELIMINATOR']
        self.parse_uniport_needed = cfg['IS_UNIPORT_NEEDS_PARSE']

    @staticmethod
    def parse_single_sample_uniport(cfg, full_uniport):
        pattern = cfg['UNIPORT_PARSE_PATTERN']
        i = 0
        while i < len(pattern):
            slices = full_uniport.split(str(pattern[i]))
            i += 1
            if int(pattern[i]) < len(slices):
                full_uniport = slices[int(pattern[i])]
            else:
                print(f"problem parsing uniport: {full_uniport}")
                return ""
            i += 1
        return full_uniport

    def parse_uniport(self, full_uniport_a, full_uniport_b, cfg):
        if self.parse_uniport_needed:
            if cfg['OUT_FILE_NAME'] == "max_linker_xl":
                return XlParser.parse_maxlinker_uniport(full_uniport_a, full_uniport_b)
            uniport_a = XlParser.parse_single_sample_uniport(cfg, full_uniport_a)
            if full_uniport_a == full_uniport_b:
                return uniport_a, uniport_a
            else:
                uniport_b = XlParser.parse_single_sample_uniport(cfg, full_uniport_b)
                return uniport_a, uniport_b
        return full_uniport_a, full_uniport_b

    @staticmethod
    def parse_maxlinker_uniport(line_a, line_b):
        uniports_a = line_a.split(';')
        uniports_b = line_b.split(';')
        uni_a_set = set(uniports_a)
        uni_b_set = set(uniports_b)
        intersect = list(uni_a_set.intersection(uni_b_set))
        if len(intersect) > 0:
            intersect = sorted(intersect, key=lambda x: x.split('-')[0])
            intersect = sorted(intersect, key=lambda x: len(x.split('-')))
            return intersect[0], intersect[0]

        uniports_a = sorted(uniports_a, key=lambda x: x.split('-')[0])
        uniports_a = sorted(uniports_a, key=lambda x: len(x.split('-')))

        uniports_b = sorted(uniports_b, key=lambda x: x.split('-')[0])
        uniports_b = sorted(uniports_b, key=lambda x: len(x.split('-')))

        return uniports_a[0], uniports_b[0]
        # sliced_uniports = [item.split('-')[0] for item in uniports]
        # i = 1
        # while i < len(sliced_uniports):
        #     if sliced_uniports[0] != sliced_uniports[i]:
        #         return ""
        #     i += 1
        # return sliced_uniports[0]

    @staticmethod
    def parse_maxlinker_res_num(line, uniport):
        candidates = line.split(';')
        for c in candidates:
            if c.split('_')[0] == uniport:
                unipot_pos = c.split('_')
                return unipot_pos[1]
        # candidates = sorted(candidates, key=lambda x: x.split('_')[0])
        # candidates = sorted(candidates, key=lambda x: len(x.split('-')))
        # candidate = candidates[0]
        # unipot_pos = candidate.split('_')
        # only_uniport = unipot_pos[0].split('-')
        # if len(only_uniport) == 1:
        # return unipot_pos[1]
        return '-1'

    @staticmethod
    def parse_res_num(cols, cfg, uniport_a, uniport_b):
        if cfg['RES_A_POS'] >= 0:
            if not cfg['IS_RES_POS_NEEDS_PARSE']:
                res_num_a = cols[cfg['RES_A_POS']]
                res_num_b = cols[cfg['RES_B_POS']]
            elif cfg['OUT_FILE_NAME'] == "max_linker_xl":
                res_num_a = XlParser.parse_maxlinker_res_num(cols[cfg['RES_A_POS']], uniport_a)
                res_num_b = XlParser.parse_maxlinker_res_num(cols[cfg['RES_B_POS']], uniport_b)
            else:
                return '-1', '-1'
        else:
            res_num_a = str(int(cols[cfg['PEP_A_POS']]) + int(cols[cfg['POS_IN_PEP_A']]) - 1)
            res_num_b = str(int(cols[cfg['PEP_B_POS']]) + int(cols[cfg['POS_IN_PEP_B']]) - 1)
        return res_num_a, res_num_b

    @staticmethod
    def parse_maxlinker_pos_in_pep(line):
        i = 0
        while i < len(line):
            if str.islower(line[i]):
                return i
            i += 1
        return -1

    @staticmethod
    def parse_pos_in_pep(line, cfg):
        if cfg['POS_IN_PEP_A'] == -1:
            return -1, -1
        if cfg['OUT_FILE_NAME'] == "max_linker_xl":
            a = XlParser.parse_maxlinker_pos_in_pep(line[cfg['PEP_A']])
            b = XlParser.parse_maxlinker_pos_in_pep(line[cfg['PEP_B']])
            return a, b
        if cfg['OUT_FILE_NAME'] == 'xiview':
            return 0, 0
        return line[cfg['POS_IN_PEP_A']], line[cfg['POS_IN_PEP_B']]

    @staticmethod
    def get_pdb_path(cfg, cols):
        if cfg['PDB_NAME'] is not None or cfg['PDB_LOC'] != -1:
            if cfg['PDB_LOC'] != -1:
                # if each line has different pdb
                pdb_path = cfg['PDB_BASE_DIR'] + cols[cfg['PDB_LOC']]
            else:
                # if all samples from same pdb
                pdb_path = cfg['PDB_BASE_DIR'] + [cfg['PDB_NAME']]
            return pdb_path
        return ""

    @staticmethod
    def get_chains(cfg, cols):
        if cfg['CHAIN_A'] != -1 and cfg['CHAIN_B'] != -1:
            return cols[cfg['CHAIN_A']], cols[cfg['CHAIN_B']]
        return "", ""

    @staticmethod
    def get_distance(cfg, cols):
        """
        If available, returns distance between the crosslinked residues
        """
        if cfg['DISTANCE'] != -1:
            return cols[cfg['DISTANCE']]
        return -1

    def parse(self, cfg, xl_file):
        xl_dict = {}
        cross_links = []
        with open(xl_file, 'r') as f:
            first_line = f.readline()
            for line in f:
                cols = np.array(line.split(self.deliminator))
                uniport_a, uniport_b = self.parse_uniport(cols[cfg['UNIPORT_A']], cols[cfg['UNIPORT_B']], cfg)
                if uniport_a == "" or uniport_b == "":
                    print(f"problem parsing uniports: {cols}")
                    continue
                res_num_a, res_num_b = XlParser.parse_res_num(cols, cfg, uniport_a, uniport_b)
                if res_num_a == -1 or res_num_b == -1:
                    print(f"problem parsing res_num: {cols}")
                    continue
                key_a = uniport_a + '+' + res_num_a
                key_b = uniport_b + '+' + res_num_b
                if key_a in xl_dict and key_b in xl_dict[key_a]:
                    print(f"duplicated xl: {key_a}")
                    continue
                pos_in_pep_a, pos_in_pep_b = XlParser.parse_pos_in_pep(cols, cfg)
                pdb_path = XlParser.get_pdb_path(cfg, cols)
                distance = XlParser.get_distance(cfg, cols)
                chain_a, chain_b = XlParser.get_chains(cfg, cols)
                if cfg['PEP_A'] != -1:
                    xl_obj = cross_link.CrossLink(cols[cfg['PEP_A']], pos_in_pep_a, res_num_a,
                                                  uniport_a, cols[cfg['PEP_B']],
                                                  pos_in_pep_b, res_num_b, uniport_b,
                                                  cfg['LINKER_TYPE'], cfg['PDB_NAME'], cfg['OUT_FILE_NAME'],
                                                  pdb_path, chain_a, chain_b, distance)
                else:
                    xl_obj = cross_link.CrossLink("", cfg['POS_IN_PEP_A'], res_num_a, uniport_a, "",
                                                  cfg['POS_IN_PEP_B'], res_num_b, uniport_b, cfg['LINKER_TYPE'],
                                                  cfg['PDB_NAME'], cfg['OUT_FILE_NAME'], pdb_path, chain_a, chain_b,
                                                  distance)
                cross_links.append(xl_obj)
                if key_a in xl_dict:
                    if key_b not in xl_dict:
                        xl_dict[key_b] = {key_a}
                    else:
                        xl_dict[key_b].add(key_a)
                    xl_dict[key_a].add(key_b)
                else:
                    xl_dict[key_a] = {key_b}
                    xl_dict[key_b] = {key_a}
        print("found cross links: ", len(cross_links))
        return cross_links

    @staticmethod
    def save_cross_links(cross_links, cfg=None, name="xl_objects_xldb_data", dir_="xl_objects/"):
        if cfg is not None:
            file_name = dir_ + "xl_objects_" + cfg['OUT_FILE_NAME']
        else:
            file_name = dir_ + name
        general_utils.save_obj(cross_links, file_name)

    @staticmethod
    def convert_old_xl_list_to_cross_link_objects():
        xl_list = general_utils.load_obj('xl_with_dup_all_cols')
        linker_types = general_utils.load_obj('linker_types')
        object_list = []
        already_created = set()
        for sample in xl_list:
            pdb_file = ''
            if sample[general_utils.XL_UNIPORT_A] == sample[general_utils.XL_UNIPORT_B]:
                if sample[general_utils.XL_PDB_TYPE_A] == general_utils.XL_AVAILABLE_IN_PDB:
                    pdb_file = sample[general_utils.XL_PDB_A]
            else:
                if sample[general_utils.XL_STRUCTURES] != '':
                    pdb_file = sample[general_utils.XL_STRUCTURES].split(',')[-1]
            if sample[general_utils.XL_DATASETS] in linker_types:
                linker = linker_types[sample[general_utils.XL_DATASETS]]
            else:
                linker = ''
            cl = cross_link.CrossLink(sample[general_utils.XL_PEP_A], sample[general_utils.XL_POS_IN_PEP_A],
                           sample[general_utils.XL_RES_NUM_A], sample[general_utils.XL_UNIPORT_A],
                           sample[general_utils.XL_PEP_B], sample[general_utils.XL_POS_IN_PEP_B],
                           sample[general_utils.XL_RES_NUM_B], sample[general_utils.XL_UNIPORT_B],
                           linker, pdb_file, sample[general_utils.XL_DATASETS])
            key1, key2 = cl.get_keys()
            if not (key1 in already_created or key2 in already_created):
                object_list.append(cl)
                already_created.add(key1)

        XlParser.save_cross_links(object_list)

    @staticmethod
    def export_xl_objects_to_csv(xl_objects, out_path):
        with open(out_path, 'w') as f:
            f.write("res_num_a,chain_a,uniport_a,pep_a,pos_in_pep_a,res_num_b,chain_b,uniport_b,pep_b,pos_in_pep_b,"
                    "distance,cb_distance,linker_type,pdb_file,origin_db,af_pae_error,pdb_path,xl_type\n")
            for obj in xl_objects:
                f.write(f"{obj.res_num_a},{obj.chain_a},{obj.uniport_a},{obj.pep_a},{obj.pos_in_pep_a},{obj.res_num_b},"
                        f"{obj.chain_b},{obj.uniport_b},{obj.pep_b},{obj.pos_in_pep_b},{obj.distance},"
                        f"{obj.cb_distance},{obj.linker_type},{obj.pdb_file},{obj.origin_db},{obj.error},"
                        f"{obj.pdb_path},{obj.xl_type},{obj.res_a_type},{obj.res_b_type}\n")

    @staticmethod
    def convert_objects_pdb_paths(new_dir, xl_objects):
        for obj in xl_objects:
            if obj.pdb_path != '':
                suff = obj.pdb_path.split('/')[-1]
                obj.pdb_path = new_dir + suff

    @staticmethod
    def import_xl_objects_from_csv(in_path, new_pdb_dir=None):
        xl_objects = []
        with open(in_path, 'r') as f:
            r = csv.reader(f, delimiter=',')
            next(r, None) # skip header
            for line in r:
                obj = cross_link.CrossLink(line[3], line[4], line[0], line[2], line[8], line[9], line[5], line[7],
                                           line[12], line[13], line[14], line[16], line[1], line[6], float(line[10]), float(line[11]),
                                           int(line[17]), float(line[15]), line[18], line[19])
                xl_objects.append(obj)
        if new_pdb_dir is not None:
            XlParser.convert_objects_pdb_paths(new_pdb_dir, xl_objects)
        return xl_objects


def run_parsing(cfg_path = '/cs/labs/dina/seanco/xl_parser/configs/parser_xl_db.yaml'):
    cfg = load_config(cfg_path)
    p = XlParser(cfg)
    xl_obj_list = p.parse(cfg, cfg['DATA_FILE'])
    return xl_obj_list
    # XlParser.save_cross_links(xl_obj_list, cfg)
    # XlParser.convert_old_xl_list_to_cross_link_objects()


# def main():
#     run_parsing()
#     # XlParser.convert_old_xl_list_to_cross_link_objects()
#
#
# if __name__ == "__main__":
#     main()
