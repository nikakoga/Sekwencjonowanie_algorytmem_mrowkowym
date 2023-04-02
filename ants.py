# pylint: disable=unsubscriptable-object
"""Moduł dla mrówków"""
from typing import Tuple
import codecs
from copy import deepcopy
from collections import defaultdict
from random import randint, choices
import numpy as np
import textdistance as td
from Original_sequence import Original_sequence, sequence_generator
from random_sequencing import regenerate_sequence, EMPTY_ELEMENT

ITERATIONS = 50
ANTS = 50
DNA_LENGTH = 500
OLIGO_LENGTH = 7
K = 10  # K best solutions from random sequencing
ITERATIONS_TO_EVAPORATE_PHEROMONES = 1
PHEROMONES_TO_EVAPORATE = 0.1
FIRST_ITERATION_MODIFIER = 3

log_file: codecs.StreamReaderWriter = None
print_enabled: bool = False


def log(text: str):
    """Wrapper na pisanie do logu"""
    if log_file is not None:
        log_file.write(text)
        log_file.write("\n")


def log_best_solution(best_solution: Tuple[float, str], original_sequence: str, iteration: int):
    """Zapisuje do pliku logu najlepsze rozwiazanie"""
    dist = td.levenshtein.distance(
        original_sequence, best_solution[1])
    log(f"Zaktualizowano najlepszy wynik po iteracji: {iteration}")
    log(f"Wartość funkcji celu: {best_solution[0]}")
    log(f"Dystans levensteina: {dist}\n")


def log_iteration(pheromone_usage_chance: float, top_k: list, iteration: int, original_sequence: str):
    """Zapisuje do pliku logu iterację"""
    log(f"\nIteracja {iteration}\n")
    log(f"Szansa na użycie feromonów {pheromone_usage_chance}")
    for i, (k, _, _, seq) in enumerate(top_k):
        dist = td.levenshtein.distance(
            original_sequence, seq)
        log(f"\t{i + 1}. Wartość funkcji celu: {k}")
        log(f"\t\tDystans Levenshteina: {dist}")


def m_print(text: str, end="\n"):
    '''Printing wrapper'''
    if print_enabled:
        print(text, end=end)


def get_target_fn(visited_times: dict, path: list) -> float:
    """obliczanie funkcji celu dla randomowych rozwiazań"""
    visited_weights = (visited_times[1] * 1 + visited_times[2] * 2 + visited_times[3] * 3) / len(path)
    # przejscia po lukach z wagą wieksza niz 1
    greater_than_1 = (visited_times[2] * 1 + visited_times[3] * 2) / len(path)
    # ignrouje kiedy takich nie bylo
    if greater_than_1 == 0:
        greater_than_1 = 1
    return visited_weights * greater_than_1


def sequencing_no_pheromones(matrix: np.array, start_row: int, original_sequence: Original_sequence,
                             seq: str, visited: set, res: defaultdict, path: list) -> Tuple[list[str], defaultdict, list]:
    '''Sekwencjonowanie bez użycia feromonów'''
    # ustawiam na 4 po to aby przejsc do czegos z wagi 3 jesli nie bedzie nic lepszego
    visited.add(start_row)
    possibilities = []
    all_neighbours = []
    cur_idx = None
    for idx, col in enumerate(matrix[start_row]):
        if col != EMPTY_ELEMENT:
            # jesli to co w kolumnie ma wage mniejsza niz aktualnie najlepsza opcja i nie jest to 0 (czyli krawedz wlasna)
            # i o ile waga jest mniejsza lub rowna niz to ile nt nam brakuje do sekwencji oryginalnej
            if idx not in visited:
                possibilities.append((idx, col))  # jesli tam nie bylam to jest to jedna z opcji
            all_neighbours.append((idx, col))  # a generalnie jest to jeden z sasiadow

    if len(possibilities) > 0:
        cur_idx, cur_sim = choices(population=possibilities, k=1)[0]  # jesli w ogole mamy jakies opcje nieodwiedzone to z nich losuje

    if cur_idx is None:
        cur_idx, cur_sim = choices(population=all_neighbours, k=1)[0]  # jesli nie mamy to losuje jakiegokolwiek sasiada
        # choices zwraca liste a ja chce tylko pojedynczy element wiec biore pierwsza rzecz z tej listy

    path.append(cur_idx)  # dodaje to co wylosowalam do sciezki
    res[cur_sim] += 1  # dodaje w zaleznosci od tego po jakiej wadze przeszlam, ze ta waga jest teraz uzyta o raz wiecej
    # dodaje wylosowana opcje
    nucleotid_to_add = original_sequence.oligo_list[cur_idx]  # dodaje aktualny wierzcholek
    seq += (nucleotid_to_add[-cur_sim:])  # a dokladnie tyle z niego od konca ile wynosi waga

    if len(seq) >= len(original_sequence.sequence):  # koncze dzialanie algorytmu gdy odtworzona sekwencja ma dlugosc wieksza lub rowna oryginalnej
        return seq, res, path
    return sequencing_no_pheromones(matrix, cur_idx, original_sequence, seq, visited, res, path)


def sequencing_with_pheromones(matrix: np.array, start_row: int, original_sequence: Original_sequence,
                               seq: str, visited: set, res: defaultdict, path: list, pheromones: np.array) -> Tuple[list[str], defaultdict, list]:
    '''Sekwencjonowanie z użyciem feromonów'''
    visited.add(start_row)
    possibilities = []
    all_neighbours = []
    cur_idx = None
    for idx, col in enumerate(matrix[start_row]):
        if col != EMPTY_ELEMENT:
            # jesli to co w kolumnie ma wage mniejsza niz aktualnie najlepsza opcja i nie jest to 0 (czyli krawedz wlasna)
            # i o ile waga jest mniejsza lub rowna niz to ile nt nam brakuje do sekwencji oryginalnej
            if idx not in visited:
                possibilities.append((idx, col))
            all_neighbours.append((idx, col))

    if len(possibilities) > 0:
        cur_idx, cur_sim = choices(population=possibilities, weights=map(lambda x: pheromones[x[0]][x[0]], possibilities), k=1)[0]

    if cur_idx is None:
        cur_idx, cur_sim = choices(population=all_neighbours, weights=map(lambda x: pheromones[x[0]][x[0]], all_neighbours), k=1)[0]

    path.append(cur_idx)
    res[cur_sim] += 1
    # dodaje wylosowany
    nucleotid_to_add = original_sequence.oligo_list[cur_idx]
    seq += (nucleotid_to_add[-cur_sim:])

    if len(seq) >= len(original_sequence.sequence):  # koncze dzialanie algorytmu gdy odtworzona sekwencja ma dlugosc wieksza lub rowna oryginalnej
        return seq, res, path
    return sequencing_with_pheromones(matrix, cur_idx, original_sequence, seq, visited, res, path, pheromones)


def regenerate_sequence_no_pheromones(dna: Original_sequence, startrow: int = 0) -> Tuple[list[str], defaultdict, list]:
    """Regenerowanie sekwencji losowo"""
    mat = np.array(dna.oligo_matrix, dtype=int)
    (seq, res, path) = sequencing_no_pheromones(
        mat, startrow, dna, dna.oligo_list[0], set(), defaultdict(lambda: 0), [startrow])  # tworzenie hashsetu
    return seq, res, path


def regenerate_sequence_with_pheromones(dna: Original_sequence, pheromones: np.array, startrow: int = 0) -> Tuple[list[str], defaultdict, list]:
    """Regenerowanie sekwencji z uzyciem feromonow"""
    mat = np.array(dna.oligo_matrix, dtype=int)
    (seq, res, path) = sequencing_with_pheromones(
        mat, startrow, dna, dna.oligo_list[0], set(), defaultdict(lambda: 0), [startrow], pheromones)  # tworzenie hashsetu
    return seq, res, path


def regenerate_sequence_with_ants(oligo_length: int = OLIGO_LENGTH,
                                  dna_length: int = DNA_LENGTH,
                                  ants_c: int = ANTS,
                                  iterations: int = ITERATIONS,
                                  positive: int = 0,
                                  negative: int = 0,
                                  evaporate_coeff: int = 1):
    """odtwarza sekwe uzywajac algorytmu mrowek"""
    increment_pheromone_usage = 100/iterations
    generated_dna = sequence_generator(oligo_length, dna_length)
    dna = Original_sequence(generated_dna, negative, positive, oligo_length, dna_length)
    dna.generate_oligo()
    dna.create_oligo_matrix()
    original_sequence_str = str.join('', dna.sequence)

    pheromones = np.zeros((len(dna.oligo_list), len(dna.oligo_list)))
    pheromone_usage_chance = 0

    best_solution: Tuple[float, str] = None
    for i in range(iterations):
        results = []
        ants_count = ants_c
        if i == 0:
            ants_count *= FIRST_ITERATION_MODIFIER
        for j in range(ants_count):
            m_print(f"\rIteracja {i + 1}/{iterations}, mrówka {j + 1}/{ants_count}", end="")
            if i == 0:  # pierwsza iteracja - algorytm losowego sekwencjonowania
                (seq, visited_times, path) = regenerate_sequence(deepcopy(dna), randomize=True)
                target_fn_val = get_target_fn(visited_times, path)
                results.append((target_fn_val, visited_times, path, str.join('', seq)))
            else:
                random_int = randint(0, 100)
                use_pheromone = random_int < pheromone_usage_chance
                if use_pheromone:
                    m_print(" - używa feromonów", end="")
                    (seq, visited_times, path) = regenerate_sequence_no_pheromones(deepcopy(dna))
                    target_fn_val = get_target_fn(visited_times, path)
                    results.append((target_fn_val, visited_times, path, str.join('', seq)))
                else:
                    m_print(" - nie używa feromonów", end="")
                    (seq, visited_times, path) = regenerate_sequence_with_pheromones(deepcopy(dna), pheromones)
                    target_fn_val = get_target_fn(visited_times, path)
                    results.append((target_fn_val, visited_times, path, str.join('', seq)))

        k_count = K
        if i == 0:
            k_count *= FIRST_ITERATION_MODIFIER

        results.sort(key=lambda x: x[0], reverse=False)  # sort by weight sum
        results = results[:k_count * 2]
        results.sort(key=lambda x: len(set(x[2]))/len(dna.oligo_matrix), reverse=True)  # sort by coverage
        results = results[:k_count]

        for k, (_, _, path, _) in enumerate(results):  # top K wyników
            for verticle in path:
                pheromones[verticle][verticle] += k_count - k
        top_k_best_target_fn_val = results[0][0]

        if best_solution is None or (best_solution is not None and best_solution[0] > top_k_best_target_fn_val):
            seq = results[0][3]
            best_solution: Tuple[float, str] = (top_k_best_target_fn_val, seq)
            log_best_solution(best_solution, original_sequence_str, i + 1)
        pheromone_usage_chance += increment_pheromone_usage
        log_iteration(pheromone_usage_chance, results, i + 1, original_sequence_str)
        if ((i + 1) % ITERATIONS_TO_EVAPORATE_PHEROMONES) == 0:
            for i, row in enumerate(pheromones):
                for j, col in enumerate(row):
                    pheromones[i][j] = col * ((1.0 - (PHEROMONES_TO_EVAPORATE * evaporate_coeff)))

    m_print("")
    log("\n========================================\n")
    log(f"Feromony na koniec: {pheromones}")
    log(f"Suma odwiedzonych wag najlepszego rozwiązania: {best_solution[0]}")
    log(f"Najlepsze rozwiązanie: {best_solution}")
    dist = td.levenshtein.distance(
        original_sequence_str, best_solution[1])
    log(f"oryginalna: {original_sequence_str}\nodtworzona: {best_solution[1]}\ndystans Levenshteina: {dist}")
    return pheromones, best_solution, original_sequence_str


def main():
    '''"Main" fn'''
    regenerate_sequence_with_ants()


if __name__ == "__main__":
    print_enabled = True
    with codecs.open("log.txt", "w", "utf-8") as log_file:
        main()
