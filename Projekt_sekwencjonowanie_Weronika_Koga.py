from collections import defaultdict
import numpy as np
import textdistance as td
from Original_sequence import Original_sequence, sequence_generator
from random_sequencing import not_perfect_sequencing


if __name__ == "__main__":
    generated_DNA = sequence_generator(7, 100)
    DNA = Original_sequence(generated_DNA, 0, 0, 7, 100)
    DNA.generate_oligo()
    DNA.create_oligo_matrix()

    mat = np.array(DNA.oligo_matrix, dtype=int)
    (reconstructed_sequence_str, _) = not_perfect_sequencing(
        mat, 0, DNA, DNA.oligo_list[0], set(), defaultdict(lambda: 0), [0])  # tworzenie hashsetu
    # połączenie elementów listy w stringa aby wyświetlić
    print(reconstructed_sequence_str)
    print(f"len(reconstruct) = {len(reconstructed_sequence_str)}")
    original_sequence_str = str.join('', DNA.sequence)
    print(original_sequence_str)
    print(f"len(DNA.sequence) = {len(DNA.sequence)}")

    similarity = td.levenshtein.similarity(
        original_sequence_str, reconstructed_sequence_str)
    print(f"DNA {100} nt\n")
    print(f"Ilosc oligo {100*7}")
    print(f"Odleglosc Levensteina {similarity}\n")
    #print(f"ZGODNOŚĆ: {similarity/len(original_sequence_str)*100}%")
