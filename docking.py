"""Docking

Script which prepares ligands and runs Autodock vina
"""
import subprocess
import os
from tqdm import tqdm


def main():
    file_name_lipinski ="results-lipinski.csv"
    pdbqt_path = "obabel_ligands_pdbqt/"
    smi_path = "obabel_ligands_smi/"
    results_path = 'docking_results/'

    protein_box = {"7nn0_A_rec.pdbqt": {"size_x": 20, "size_y": 24, "size_z": 25,
                                        "center_x":110, "center_y": -5, "center_z": 48},
                    "6zsl_A_rec.pdbqt": {"size_x": 20, "size_y": 20, "size_z": 20,
                                        "center_x":-15, "center_y": 10, "center_z": -80},
                    "5rme_A_rec.pdbqt": {"size_x": 20, "size_y": 25, "size_z": 25,
                                        "center_x":-12, "center_y": 15, "center_z": -78},
                    "5rm2_A_rec.pdbqt": {"size_x": 20, "size_y": 20, "size_z": 20,
                                        "center_x":-10, "center_y": 10, "center_z": -80},
                    "7nio_A_rec.pdbqt": {"size_x": 20, "size_y": 20, "size_z": 20,
                                        "center_x":-33, "center_y": 14, "center_z": -43}}

    make_smi_files(smi_path, file_name_lipinski)
    make_pdbqt_files(pdbqt_path, smi_path)
    run_vina(results_path, protein_box, pdbqt_path)


def make_smi_files(smi_dir:str, file_name:str):
    """Function which for every ligand creates separate .smi file

    :param smi_dir: Path to which smi files will be saved
    :type smi_dir: str
    :param file_name: File with ligands which follow the Lipinski rule
    :type file_name: str
    """
    if not os.path.exists(smi_dir):
        os.makedirs(smi_dir)

    with open(file_name, "r", encoding="UTF-8") as f_in:
        lines = f_in.readlines()
        for line in lines[1:502]:
            words = line.split(",")
            if not os.path.isfile(smi_dir + "/" + words[-2]):
                with open(smi_dir + words[-2]+".smi", "w", encoding="UTF-8") as f_out:
                    f_out.write(words[1])
                with open("names.txt", "a") as f_names:
                    f_names.write(words[-2] + ".pdbqt" + "\n" )


def make_pdbqt_files(pdbqt_dir:str, smi_dir:str):
    """Function which for every ligand saved in smi file creates pdbqt file

    :param pdbqt_dir: Path to which pdbqt files will be saved
    :type pdbqt_dir: str
    :param smi_dir: Path to folder with .smi files
    :type smi_dir: str
    """
    if not os.path.exists(pdbqt_dir):
        os.makedirs(pdbqt_dir)

    for file in os.listdir(smi_dir):
        command = "obabel "+ smi_dir + file + " -O " + pdbqt_dir + file[:-3] + "pdbqt --gen3d"
        subprocess.run(command, shell = True)



def run_vina(path_results:str, prot_box:dict, pdbqt_dir:str):
    """Function which runs Autodock vina for each of
    the 5 proteins and every ligand

    :param path_results: Path to folder which stores results
    :type path_results: str
    :param prot_box: Dictionary with boz sizes and placement for every protein
    :type prot_box: dict
    :param pdbqt_dir: Path to folder which stores ligands
    :type pdbqt_dir: str
    """
    if not os.path.exists(path_results):
        os.makedirs(path_results)

    for protein in prot_box.keys():
        print("Docking ligands for protein: ", protein[:-6])
        if not os.path.exists(path_results+protein[:-6]):
            os.makedirs(path_results+protein[:-6])
        for ligand in tqdm(os.listdir(pdbqt_dir)):
            if ligand.endswith(".pdbqt"):
                if not os.path.exists(protein[:-6] + "/" + ligand[:-6]):
                    command = "vina --receptor {} --ligand {} --size_x {} --size_y {} " \
                        "--size_z {} --center_x {} --center_y {} --center_z {} --seed 123 " \
                            "--log {}.txt --out {}_out.pdbqt ".format("proteins/" + protein,
                                                                     pdbqt_dir+ligand,
                                                                     prot_box[protein]["size_x"],
                                                                     prot_box[protein]["size_y"],
                                                                     prot_box[protein]["size_z"],
                                                                     prot_box[protein]["center_x"],
                                                                     prot_box[protein]["center_y"],
                                                                     prot_box[protein]["center_z"],
                                                                     path_results + protein[:-6] + "/" + ligand[:-6],
                                                                     path_results + protein[:-6] + "/" + ligand[:-6])
                    subprocess.run(command, shell=True)


if __name__ == "__main__":
    main()
