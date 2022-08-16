import warnings
from Bio.PDB.Atom import Atom
from enum import Enum
from os import environ
from numpy import array


protein_letters = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS",
    "L": "LEU", "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL",
    "W": "TRP", "Y": "TYR",
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F", "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K",
    "LEU": "L", "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R", "SER": "S", "THR": "T", "VAL": "V",
    "TRP": "W", "TYR": "Y"
}


def get_amino_acid(amino_acid):
    """
    Get correct amino acid data from segment.

    :param amino_acid: original amino acid name (maybe wrong or other char type).
    :type amino_acid: str

    :return: right amino acid name.
    :rtype: str
    """
    amino_acid = amino_acid.upper()
    if amino_acid in protein_letters:
        return protein_letters[amino_acid]
    else:
        return "-"


def get_peptides(protein_path, minimum_length=1):
    """
    Obtain peptides from the protein file.

    :param protein_path: path for loading protein file.
    :type protein_path: str

    :param minimum_length: minimum length of peptide.
    :type minimum_length: int

    :returns: peptide chains, their short data, and additional residue data if the file has.
    :rtype: dict or None

    .. note::
        Reference [1]: Andrei Kouranov, Lei Xie, Joanna de la Cruz, and others (2006) Nucleic Acids Res.
        Reference [2]: John Jumper, Richard Evans, Alexander Pritzel, and others (2021). Nature.
    """
    total_peptides, modified_residues, available_peptides, positions = dict(), dict(), dict(), dict()
    prepare_finish = False
    with open(protein_path, "r") as pdb_file:
        current_residue_index, first_atom_index = "", -1
        for line in pdb_file.readlines():
            # find the resolution in the PDB file.
            if "REMARK   2 RESOLUTION." in line:
                resolution = float(line[23:line.index("ANGSTROMS") - 2])

            # find all the residues.
            header = line[:line.index(" ")]
            if header == "SEQRES":
                chain_index = line[11]
                residues = list(map(get_amino_acid, filter(lambda r: len(r) == 3, line[19:70].split(" "))))
                if chain_index in total_peptides:
                    total_peptides[chain_index] += "".join(residues)
                else:
                    total_peptides[chain_index] = "".join(residues)

            # find the modification in residues if the protein has.
            elif header == "MODRES":
                base_residue = get_amino_acid(line[24:27])
                if base_residue != "-":
                    modified_residues[(line[16], line[12:15], int(line[18:22]) - 1)] = base_residue

            # find N-CA-C skeleton (protein backbone) in atom data.
            elif header == "ATOM" or header == "HETATM":
                # finish the peptide data first.
                if not prepare_finish:
                    if len(modified_residues) > 0:
                        for (chain_index, changed_residue, position), base_residue in modified_residues.items():
                            if position < len(total_peptides[chain_index]):
                                sequence = list(total_peptides[chain_index])
                                sequence[position] = base_residue
                                total_peptides[chain_index] = "".join(sequence)

                    temp_peptides, temp_modified_residues = {}, {}
                    for chain_index, peptide in total_peptides.items():
                        if ("-" not in peptide) and (len(peptide) > 0):
                            temp_peptides[chain_index] = peptide
                            for information, base_residue in modified_residues.items():
                                if information[0] == chain_index:
                                    temp_modified_residues[information] = base_residue
                    total_peptides, modified_residues = temp_peptides, temp_modified_residues
                    prepare_finish = True

                # read atom data from file.
                atom, atom_index, residue = line[12:16].strip(), int(line[6:11]), get_amino_acid(line[17:20])
                chain_index, residue_index, insert_flag = line[21], int(line[22:26]), line[26]
                position = [float(line[30:38]), float(line[38:46]), float(line[46:54])]

                # ignore invalid chain.
                if chain_index not in list(total_peptides.keys()):
                    continue

                # check residue from peptide data.
                for (changed_chain_index, _, changed_position), base_residue in modified_residues.items():
                    if chain_index == changed_chain_index and residue_index == changed_position:
                        residue = base_residue
                        break

                # add new chain index.
                if chain_index not in positions:
                    positions[chain_index] = []
                    # starts[chain_index] = residue_index
                    available_peptides[chain_index] = ""

                # check atom need to be collected.
                if current_residue_index != str(residue_index) + insert_flag:
                    current_residue_index = str(residue_index) + insert_flag
                    first_atom_index = atom_index
                    available_peptides[chain_index] += residue
                if atom_index - first_atom_index <= 2:  # N-CA-C skeleton
                    positions[chain_index].append(position)

            # find end flag if reading all the available atoms.
            elif header == "TER" and total_peptides.keys() == positions.keys():
                break

    for chain_index in positions.keys():
        positions[chain_index] = array(positions[chain_index])

    peptide_set, position_set = {}, {}
    for chain_index in total_peptides.keys():
        if chain_index in available_peptides:
            if available_peptides[chain_index] == total_peptides[chain_index]:
                name = chain_index + "-0"
                peptide_set[name], position_set[name] = total_peptides[chain_index], positions[chain_index]
            else:
                location_a, temp_peptide = 0, ""
                for location_t in range(len(total_peptides[chain_index])):
                    if (location_a < len(available_peptides[chain_index])) and \
                            available_peptides[chain_index][location_a] == total_peptides[chain_index][location_t]:
                        temp_peptide += available_peptides[chain_index][location_a]
                        location_a += 1
                    else:
                        temp_peptide += "-"

                # Ignore this peptide because it catches many mistakes.
                if temp_peptide.count("-") > len(total_peptides[chain_index]) * 0.1:
                    continue

                peptide_parts = list(filter(lambda p: len(p) > minimum_length, temp_peptide.split("-")))
                for part_index, peptide_part in enumerate(peptide_parts):
                    start_index = total_peptides[chain_index].index(peptide_part)
                    last_length = sum([len(p) for p in peptide_part[: part_index]]) if part_index > 1 else 0
                    name = chain_index + "-" + str(start_index)
                    peptide_set[name] = peptide_part
                    position_set[name] = positions[chain_index][last_length * 3: (last_length + len(peptide_part)) * 3]

    return {"total": total_peptides, "modified": modified_residues, "usable": peptide_set,
            "positions": position_set, "resolution": resolution}


def get_structures(protein_path, peptides, modified_residues):
    """
    Get secondary structure data of protein file.

    :param protein_path: path for loading protein file.
    :type protein_path: str

    :param peptides: obtained peptide information in the corresponding protein path.
    :type peptides: dict

    :param modified_residues: modified residues in the corresponding protein path.
    :type modified_residues:dict

    :returns: secondary structure data.
    :rtype dict

    .. note::
        Reference: Andrei Kouranov, Lei Xie, Joanna de la Cruz, and others (2006) Nucleic Acids Res.
    """
    secondaries = {"HELIX": [], "SHEET": [], "TURN": [], "SSBOND": [],
                   "LINK": [], "HYDBND": [], "SLTBRG": [], "CISPEP": []}

    with open(protein_path, "r") as pdb_file:
        for line in pdb_file.readlines():
            header = line[:line.index(" ")]
            found = False
            if header == "HELIX":
                chain_index_1, chain_index_2 = line[19], line[31]
                location_index_1, location_index_2 = int(line[21:25]), int(line[33:37])
                residue_1, residue_2 = get_amino_acid(line[15:18]), get_amino_acid(line[27:30])
                found = True
            elif header == "SHEET":
                chain_index_1, chain_index_2 = line[21], line[32]
                location_index_1, location_index_2 = int(line[22:26]), int(line[33:37])
                residue_1, residue_2 = get_amino_acid(line[17:20]), get_amino_acid(line[28:31])
                found = True
            elif header == "TURN":
                chain_index_1, chain_index_2 = line[19], line[30]
                location_index_1, location_index_2 = int(line[21:24]), int(line[31:35])
                residue_1, residue_2 = get_amino_acid(line[15:18]), get_amino_acid(line[26:29])
                found = True
            elif header == "SSBOND":
                chain_index_1, chain_index_2 = line[15], line[29]
                location_index_1, location_index_2 = int(line[17:21]), int(line[31:35])
                residue_1, residue_2 = get_amino_acid(line[11:14]), get_amino_acid(line[25:28])
                found = True
            elif header == "LINK":
                chain_index_1, chain_index_2 = line[21], line[51]
                location_index_1, location_index_2 = int(line[22:26]), int(line[52:56])
                residue_1, residue_2 = get_amino_acid(line[17:20]), get_amino_acid(line[47:50])
                found = True
            elif header == "HYDBND":
                chain_index_1, chain_index_2 = line[21], line[52]
                location_index_1, location_index_2 = int(line[22:27]), int(line[53:58])
                residue_1, residue_2 = get_amino_acid(line[17:20]), get_amino_acid(line[48:51])
                found = True
            elif header == "SLTBRG":
                chain_index_1, chain_index_2 = line[21], line[51]
                location_index_1, location_index_2 = int(line[22:26]), int(line[52:56])
                residue_1, residue_2 = get_amino_acid(line[17:20]), get_amino_acid(line[47:50])
                found = True
            elif header == "CISPEP":
                chain_index_1, chain_index_2 = line[15], line[29]
                location_index_1, location_index_2 = int(line[17:21]), int(line[31:35])
                residue_1, residue_2 = get_amino_acid(line[11:14]), get_amino_acid(line[25:28])
                found = True

            if found:
                # check residue from peptide data.
                for (changed_chain_index, _, changed_position), base_residue in modified_residues.items():
                    if chain_index_1 == changed_chain_index and location_index_1 == changed_position:
                        residue_1 = base_residue
                        break
                for (changed_chain_index, _, changed_position), base_residue in modified_residues.items():
                    if chain_index_2 == changed_chain_index and location_index_2 == changed_position:
                        residue_2 = base_residue
                        break

                if residue_1 == "-" or residue_2 == "-":
                    continue

                check_1, check_2 = False, False
                for bias in range(3):
                    if not check_1:
                        try:
                            if peptides[chain_index_1][location_index_1 - bias] == residue_1:
                                location_index_1 -= bias
                                check_1 = True
                        except IndexError:
                            pass
                        except KeyError:
                            pass

                        try:
                            if peptides[chain_index_1][location_index_1 + bias] == residue_1:
                                location_index_1 += bias
                                check_1 = True
                        except IndexError:
                            pass
                        except KeyError:
                            pass

                    if not check_2:
                        try:
                            if peptides[chain_index_2][location_index_2 - bias] == residue_2:
                                location_index_2 -= bias
                                check_2 = True
                        except IndexError:
                            pass
                        except KeyError:
                            pass

                        try:
                            if peptides[chain_index_2][location_index_2 + bias] == residue_2:
                                location_index_2 += bias
                                check_2 = True
                        except IndexError:
                            pass
                        except KeyError:
                            pass

                if check_1 and check_2:
                    secondaries[header].append((chain_index_1, location_index_1, chain_index_2, location_index_2))

            if header == "ATOM":
                break

    return secondaries


# Type to save token PDB file.
# 'average' is to save the average individual calculated from all the individual.
# 'individual' is to save files of each individual.
class MotifFileType(Enum):
    average = 0
    individual = 1


def set_motif_file(motif_segment, motif_structures, cluster_index=None,
                   token_file_type=MotifFileType.average, save_path=None):
    """
    Set the motif file.

    :param motif_segment: token segment.
    :type motif_segment: str

    :param motif_structures: token structure groups.
    :type motif_structures: list

    :param cluster_index: cluster index of a token.
    :type cluster_index: int

    :param token_file_type: token type, average of cluster or all in cluster.
    :type token_file_type: pt.converter.TokenFileType

    :param save_path: path to save protein file.
    :type save_path: str or None

    :return: file path(s) of protein.
    :rtype: str or list
    """
    warnings.filterwarnings("ignore")
    if save_path is None:
        save_path = environ["TMP"].replace("\\", "/") + "/"
    form = ["N", "CA", "C"]
    if token_file_type == MotifFileType.average:
        average_position_data = motif_structures[0]
        for position_data in motif_structures[1:]:
            average_position_data += position_data
        average_position_data = average_position_data / len(motif_structures)
        file_path = save_path + str(len(motif_segment)) + "_" + motif_segment + "_A_"
        file_path += ("C_" + str(cluster_index + 1) if cluster_index else "") + ".pdb"
        with open(file_path, "w") as pdb_file:
            for position_index, position in enumerate(average_position_data):
                atom = Atom(name=form[position_index % len(form)], serial_number=position_index + 1,
                            fullname=" " + form[position_index % len(form)].ljust(3), altloc=" ",
                            bfactor=0.0, coord=position.tolist(), occupancy=1.0)
                residue = protein_letters[motif_segment[position_index // len(form)]].upper()
                line = write_atom_line(atom=atom, atom_number=position_index + 1, residue=residue, chain_index="A")
                pdb_file.write(line)
            pdb_file.write("TER\nEND")
        return file_path
    elif token_file_type == MotifFileType.individual:
        file_paths = []
        for file_index, position_data in enumerate(motif_structures):
            file_path = save_path + str(len(motif_segment)) + "_" + motif_segment + "_I_"
            file_path += ("C_" + str(cluster_index + 1) + "_") if cluster_index else ""
            file_path += str(file_index + 1) + ".pdb"
            with open(file_path, "w") as pdb_file:
                for position_index, position in enumerate(position_data):
                    atom = Atom(name=form[position_index % len(form)], serial_number=position_index + 1,
                                fullname=" " + form[position_index % len(form)].ljust(3), altloc=" ",
                                bfactor=0.0, coord=position.tolist(), occupancy=1.0)
                    residue = protein_letters[motif_segment[position_index // len(form)]].upper()
                    line = write_atom_line(atom=atom, atom_number=position_index + 1, residue=residue, chain_index="A")
                    pdb_file.write(line)
                pdb_file.write("TER\nEND")
            file_paths.append(file_path)
        return file_paths
    else:
        raise ValueError("No such type of token file!")


def write_atom_line(atom, atom_number, residue, chain_index):
    """
    Write the atom line as PDB file.

    :param atom: atom data.
    :type atom: Bio.PDB.Atom.Atom

    :param atom_number: number of atom.
    :type atom_number: int

    :param residue: residue name.
    :type residue: str

    :param chain_index: chain index in the protein.
    :type chain_index: str

    :return: the atom line.
    :rtype str
    """
    record_type = "ATOM  "
    element = atom.element.strip().upper().rjust(2)
    x, y, z = atom.get_coord()
    line = "%s%5i %-4s%c%3s %c%4i%c   %8.3f%8.3f%8.3f%s%6.2f      %4s%2s%2s\n" % \
           (record_type, atom_number, atom.get_fullname(), " ", residue, chain_index, (atom_number - 1) // 3, " ",
            x, y, z, "  1.00", 0.0, " ", element, "  ")
    return line
