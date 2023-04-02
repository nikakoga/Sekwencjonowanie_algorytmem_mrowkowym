import random
import math

EMPTY_ELEMENT = -1


def sequence_generator(one_oligo_length, size):
    generated_DNA = [random.choice('ATCG')for _ in range(size)]
    return generated_DNA


class Original_sequence:

    def __init__(self, sequence, negative, positive, one_oligo_length, size):
        self.sequence = sequence
        self.negative = negative
        self.positive = positive
        self.one_oligo_length = one_oligo_length
        self.oligo_list = []
        self.oligo_matrix = [[]]
        self.size = size

    def generate_oligo(self):
        for i, x in enumerate(self.sequence):
            odczyt = self.sequence[i:i+self.one_oligo_length]
            if odczyt not in self.oligo_list:
                self.oligo_list.append(odczyt)

    def show_oligo(self):
        print(self.oligo_list)
        liczba_oligo = len(self.oligo_list)
        print(f"Ilosc oligo {liczba_oligo}")
        print(f"Pol w macierzy {liczba_oligo*liczba_oligo}")

    def add_negative_mistakes(self):
        for _ in range(self.negative):  # robie to tyle razy ile ma byc bledow negatywnych
            # losuje indeks oligonukleotydu ktory wywale z tablicy
            oligo_do_usuniecia = random.randint(0, len(self.oligo_list)-1)
            self.oligo_list.remove(oligo_do_usuniecia)  # i go wywalam

    def add_positive_mistakes(self):

        for i in range(self.positive):

            if (len(self.oligo_list) == math.pow(4, self.one_oligo_length)):
                print("Nie dodam bledow pozytywnych bo wszystkie mozliwe opcje istnieja juz w liscie oligonukleotydow\n")
                return
            sekwencja_do_dodania = [random.choice('ATCG') for _ in range(self.one_oligo_length)]

            # losuje tak dlugo sekwencje oligo az nie bedzie jej w obecnie znanych oligo
            while sekwencja_do_dodania in self.oligo_list:
                sekwencja_do_dodania = [random.choice('ATCG') for _ in range(7)]
            self.oligo_list.append(sekwencja_do_dodania)

    def create_oligo_matrix(self):
        self.oligo_matrix = [[EMPTY_ELEMENT for _ in range(
            len(self.oligo_list))] for _ in range(len(self.oligo_list))]
        licznik_1 = 0
        licznik_2 = 0
        licznik_3 = 0
        for i, row in enumerate(self.oligo_list):
            for j, column in enumerate(self.oligo_list):
                if (row == column):
                    self.oligo_matrix[i][j] = 0
                elif (row[1:] == column[:-1]):
                    self.oligo_matrix[i][j] = 1
                    licznik_1 += 1
                elif (row[2:] == column[:-2]):
                    self.oligo_matrix[i][j] = 2
                    licznik_2 += 1
                elif (row[3:] == column[:-3]):
                    self.oligo_matrix[i][j] = 3
                    licznik_3 += 1
        number_of_oligo = len(self.oligo_list)
        print(f"Ilosc oligo {number_of_oligo}")
        fields_in_matrix = number_of_oligo*number_of_oligo
        print(f"Pol w macierzy {fields_in_matrix}")
