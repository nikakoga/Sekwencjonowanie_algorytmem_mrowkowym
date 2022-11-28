import random

BLEDY_NEGATYWNE = 6
BLEDY_POZYTYWNE = 0

oryginalna_sekwencja_DNA = [random.choice('ATCG') for _ in range(210)]
odczyt = ""
oligo = []
liczba_oligo = 0

for i, x in enumerate(oryginalna_sekwencja_DNA):
    print(x,i)
    odczyt = odczyt + x
    if ((i+1)%7 == 0):
        if odczyt not in oligo:
            oligo.append(odczyt)
        odczyt = ""

print(oligo)
print(len(oligo))

def tworzenie_bledow_negatywnych(oligo):
    for i in range(BLEDY_NEGATYWNE):
        random.randint(0,len(oligo)-1)
