'''Random sequencing functions'''
from collections import defaultdict
from typing import Callable, Tuple
from random import randint
import numpy as np
from Original_sequence import Original_sequence, EMPTY_ELEMENT

MIN_SIMILAR = 3
MAX_SIMILAR = 0
MAX_ITER = 50
OLIGO_LEN = 8


def not_perfect_sequencing(matrix: np.array, start_row: int, original_sequence: Original_sequence,
                           seq: str, visited: set, res: defaultdict, path: list, randomize=True) -> Tuple[str, defaultdict, list]:
    '''Generates random sequence'''
    cur_sim = MIN_SIMILAR + 1
    # ustawiam na 4 po to aby przejsc do czegos z wagi 3 jesli nie bedzie nic lepszego
    ori_seq_len = len(original_sequence.sequence)
    visited.add(start_row)
    cur_idx = None
    if randomize:
        possibilities = []
        for idx, col in enumerate(matrix[start_row]):
            if col != EMPTY_ELEMENT:
                # jesli to co w kolumnie ma wage mniejsza niz aktualnie najlepsza opcja i nie jest to 0 (czyli krawedz wlasna)
                # i o ile waga podobienstwa jest mniejsza lub rowna niz to ile nt nam brakuje do sekwencji oryginalnej
                if (cur_sim > col > MAX_SIMILAR and col <= ori_seq_len - len(seq) and idx not in visited):
                    cur_sim = col
                    possibilities = [idx]
                elif cur_sim == col and col <= ori_seq_len - len(seq) and idx not in visited:
                    possibilities.append(idx)
        if len(possibilities) > 0:
            random_possibility = randint(0, len(possibilities) - 1)
            cur_idx = possibilities[random_possibility]
    else:
        for idx, col in enumerate(matrix[start_row]):
            if col != EMPTY_ELEMENT:
                # jesli to co w kolumnie ma wage mniejsza niz aktualnie najlepsza opcja i nie jest to 0 (czyli krawedz wlasna)
                # i o ile waga podobienstwa jest mniejsza lub rowna niz to ile nt nam brakuje do sekwencji oryginalnej
                if (cur_sim > col > MAX_SIMILAR and col <= ori_seq_len - len(seq) and idx not in visited):
                    cur_sim = col
                    cur_idx = idx

    if cur_idx is None:
        # sytuacja gdy nie ma sasiadow ktorych mozemy dodac jako nieodwiedzonych
        # odpalam pomocniczy algo ktory sprawdzi gdzie isc aby w koncu bylo cos nieodwiedzonego
        dfs_path = searching_for_unvisited(matrix, start_row, visited, ori_seq_len - len(seq), randomize)
        if len(dfs_path) == 0:  # sytuacja gdy ten algo nic nie znalazl bo wszedzie bylam
            return seq, res, path
        for cur_idx, cur_sim in dfs_path:
            path.append(cur_idx)
            res[cur_sim] += 1
            visited.add(cur_idx)  # oznaczam ze jest odwiedzony
            # szukam sekwencji tego oligo
            nucleotid_to_add = original_sequence.oligo_list[cur_idx]
            # w zaleznosci od wagi dodaje tyle nukleotydow od konca
            seq += (nucleotid_to_add[-cur_sim:])
            if len(seq) >= ori_seq_len:  # koncze dzialanie algorytmu gdy odtworzona sekwencja ma dlugosc wieksza lub rowna oryginalnej
                return seq, res, path
    else:
        path.append(cur_idx)
        res[cur_sim] += 1
        # dodaje najlepsza znaleziona opcje
        nucleotid_to_add = original_sequence.oligo_list[cur_idx]
        seq += (nucleotid_to_add[-cur_sim:])

    if len(seq) >= ori_seq_len:  # koncze dzialanie algorytmu gdy odtworzona sekwencja ma dlugosc wieksza lub rowna oryginalnej
        return seq, res, path
    return not_perfect_sequencing(matrix, cur_idx, original_sequence, seq, visited, res, path, randomize)


def dfs(adj_mat: np.array, start: int, dfs_visited: list, stack: list,
        action: Callable[[int, list], None], randomize: bool,
        max_iter: int = None, iter_count: int = 0):
    '''Depth-First Search on adjacency matrix'''  # Callable jest do typowania
    dfs_visited[start] = True  # miejsce od ktorego wywoluje dfs ustawiam jako odwiedzone
    # dla wszystkich kolumn z zadanego rzedu

    if randomize:
        possibilities = []
        for i in range(len(adj_mat[start])):
            adj_element = adj_mat[start][i]  # dokladne polozenie w macierzy
            if adj_element != EMPTY_ELEMENT and not dfs_visited[i]:
                possibilities.append(i)

        while len(possibilities) > 0:
            random_idx = randint(0, len(possibilities)-1)
            i = possibilities[random_idx]
            possibilities.remove(i)
            # dodaje na stack nr columny oraz dokladne polozenie w macierzy
            stack.append((i, int(adj_element)))
            # dodalam na stack, czyli tam wejde wiec uznaje za odwiedzone zeby znowu nie dodac tego na stack i nie odwiedzac w kolko
            dfs_visited[i] = True
            # to co jest tą akcją (on_next_element_dfs) widać w funkcji search_for_unvisited gdzie wywoluje dfs.
            # Dzięki temu w każdym nowym elemencie ścieżki dodaje "wydajność" czyli ile jest nieodwiedzonych i jakiej wagi one są
            action(i, stack)  # licze "wydajnosc" tego elementu sciezki.
            if iter_count < max_iter or max_iter is None:
                dfs(adj_mat, i, dfs_visited, stack,
                    action, max_iter, iter_count + 1)
            stack.pop()
    else:
        # dla wszystkich kolumn z zadanego rzedu
        for i in range(len(adj_mat[start])):
            adj_element = adj_mat[start][i]  # dokladne polozenie w macierzy
            if adj_element != EMPTY_ELEMENT and not dfs_visited[i]:
                # dodaje na stack nr columny oraz dokladne polozenie w macierzy
                stack.append((i, int(adj_element)))
                # dodalam na stack, czyli tam wejde wiec uznaje za odwiedzone zeby znowu nie dodac tego na stack i nie odwiedzac w kolko
                dfs_visited[i] = True
                # to co jest tą akcją (on_next_element_dfs) widać w funkcji search_for_unvisited gdzie wywoluje dfs.
                # Dzięki temu w każdym nowym elemencie ścieżki dodaje "wydajność" czyli ile jest nieodwiedzonych i jakiej wagi one są
                action(i, stack)  # licze "wydajnosc" tego elementu sciezki.
                if iter_count < max_iter or max_iter is None:
                    dfs(adj_mat, i, dfs_visited, stack,
                        action, max_iter, iter_count + 1)
                stack.pop()


def calculate_coefficient(path: list, visited: set) -> int:
    """Calculates node stack coefficient"""
    if len(path) == 0:  # aby nie dzielic przez 0
        return 0.0
    # pomocnicza funkcja ktora pozwoli wiliczyc jak dobrym rozwiazaniem jest znaleziona sciezka

    def summation(path_element):
        # przypisuje wartosc krotki do 2 osobnych zmiennych. Pierwsze z krotki to indeks, drugie podobienstwo
        col, col_similarity = path_element
        if col not in visited:
            return MIN_SIMILAR - col_similarity + 1
        return 0
    # sprawdzam jaki wspolczynnik uzyska cala sciezka. Dzieki mapie dla kazdego elementu ze sciezki wywoluje te funkcje
    similarity_sum = sum(map(summation, path))
    return similarity_sum/len(path)


def searching_for_unvisited(adj_mat: np.array, start_row: int, visited: set, length_to_add: int, randomize: bool) -> list:
    '''Searches for unvisited nodes'''
    options = []

    def on_next_element_dfs(row: int, path: list):
        # w lambdzie sprawdzam czy suma wszystkich podobieństw nieodwiedzonych nt
        # jest mniejsza lub rowna temu co jeszce mamy dodac z dfs aby nie dodać za dużo
        # jesli jest mniejsza to ok a jesli jest wieksza to nie zapisuje tej sciezki w ogole jako opcje do wyboru
        if row not in visited and sum(map(lambda x: x[1] if x not in visited else 0, path)) <= length_to_add:
            # potem ta sciezka zostanie nadpisana inna sciezka dlatego do listy opcji dodawana jest ta wlasnie znaleziona sciezka razem z jej "wydajnoscia"
            # aby potem wybrac taka z najlepszym wspolczynnikiem wydajnosci
            # to jest krotka
            options.append((calculate_coefficient(path, visited), list(path)))

    dfs(adj_mat, start_row, [False]*len(adj_mat[start_row]),  # tablica z samymi False na początku bo uruchamiam dfs pierwszy raz i w ramach tej funkcji jeszcze nic nie odwiedziłam
        [], on_next_element_dfs, randomize, max_iter=MAX_ITER)  # pusta lista to stack, akcja do wykonania to on_next_element_dfs ktory zapisze wszystko przed nastepnym dfs

    # sortuje opcje tak aby najlepsza byla na poczatku, reverse zeby sortowac malejaco tak aby najlepsza byla z przodu
    options.sort(key=lambda x: x[0], reverse=True)
    if len(options) > 0:  # jesli sa jakiekolwiek opcje
        # Biore pierwsza krotke i zwracam jej drugi element czyli najlepsza sciezke, bez wydajnosci
        return options[0][1]
    return []


def regenerate_sequence(dna_sequence: Original_sequence, startrow=0, randomize=False) -> Tuple[str, defaultdict, list]:
    '''Regenerates sequence'''
    mat = np.array(dna_sequence.oligo_matrix, dtype=int)
    (reconstructed, res, path) = not_perfect_sequencing(
        mat, startrow, dna_sequence, dna_sequence.oligo_list[0], set(), defaultdict(lambda: 0), [startrow], randomize)  # tworzenie hashsetu
    return reconstructed, res, path
