
import textdistance as td
from Original_sequence import Original_sequence


def get_distance(reconstructed_sequence_str: str, dna_sequence: Original_sequence, print_enabled=False) -> float:
    '''Gets Levenshtien distance for sequences'''
    if print_enabled:
        print(f"Zrekonstruowana sekwencja: {reconstructed_sequence_str}\n")

    original_sequence_str = str.join('', dna_sequence.sequence)

    if print_enabled:
        print(f"Oryginalna sekwencja: {original_sequence_str}\n")
        print(f"Długość zrekonstruowanej sekwencji: {len(reconstructed_sequence_str)}")
        print(f"Długość oryginalnej sekwencji: {len(original_sequence_str)}\n")
        print(f"Czy takie same: {reconstructed_sequence_str == original_sequence_str}")

    distance = td.levenshtein.distance(
        original_sequence_str, reconstructed_sequence_str)

    if print_enabled:
        print(f"DNA {dna_sequence.size} nt\n")
        print(f"Ilość oligo {dna_sequence.size * dna_sequence.one_OLIGO_LENgth}")
        print(f"Odległość Levensteina: {distance}\n")
    return distance
