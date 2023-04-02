'''Module for generating plots'''
import math
import matplotlib.pyplot as plt
from ants import regenerate_sequence_with_ants
import textdistance as td


def dna_size_plot():
    """Generates plot by dna size"""
    x = list(range(100, 800, 100))
    y = []
    for l in range(100, 800, 100):
        dist = 0
        for _ in range(3):
            _, best_solution, original_seq = regenerate_sequence_with_ants(7, l, 10, 10)
            dist += td.levenshtein.distance(best_solution[1], original_seq)
        y.append(dist/3)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(x, y)
    annotate_move_y = math.floor(max(y)/20)
    for i, txt in enumerate(y):
        ax.annotate(txt, (x[i] - 20, txt + annotate_move_y))
    plt.grid(True)
    ax.set_ylim([0, ax.get_ylim()[1] + annotate_move_y * 2])
    ax.set_ylabel("Odległość Levensteina")
    ax.set_xlabel("Długość DNA")
    plt.savefig("plots/1.png")


def iteration_count_plot():
    """Generates plot by dna size"""
    x = [10, 50, 100, 200]
    y = []
    for l in x:
        dist = 0
        for _ in range(3):
            _, best_solution, original_seq = regenerate_sequence_with_ants(7, 200, 10, l)
            dist += td.levenshtein.distance(best_solution[1], original_seq)
        y.append(dist/3)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(x, y)
    annotate_move_y = math.floor(max(y)/20)
    for i, txt in enumerate(y):
        ax.annotate(txt, (x[i] - 20, txt + annotate_move_y))
    plt.grid(True)
    ax.set_ylim([0, ax.get_ylim()[1] + annotate_move_y * 2])
    ax.set_ylabel("Odległość Levensteina")
    ax.set_xlabel("Ilość iteracji")
    plt.savefig("plots/2.png")


def positive_plot():
    """Generates plot by dna size"""
    x = [2, 5, 10, 20]
    y = []
    for l in x:
        dist = 0
        for _ in range(3):
            _, best_solution, original_seq = regenerate_sequence_with_ants(7, 200, 10, 10, positive=l)
            dist += td.levenshtein.distance(best_solution[1], original_seq)
        y.append(dist/3)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(x, y)
    annotate_move_y = math.floor(max(y)/20)
    for i, txt in enumerate(y):
        ax.annotate(txt, (x[i] - 20, txt + annotate_move_y))
    plt.grid(True)
    ax.set_ylim([0, ax.get_ylim()[1] + annotate_move_y * 2])
    ax.set_ylabel("Odległość Levensteina")
    ax.set_xlabel("Ilość błędów pozytywnych[%]")
    plt.savefig("plots/3.png")


def negative_plot():
    """Generates plot by dna size"""
    x = [2, 5, 10, 20]
    y = []
    for l in x:
        dist = 0
        for _ in range(3):
            _, best_solution, original_seq = regenerate_sequence_with_ants(7, 200, 10, 10, negative=l)
            dist += td.levenshtein.distance(best_solution[1], original_seq)
        y.append(dist/3)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(x, y)
    annotate_move_y = math.floor(max(y)/20)
    for i, txt in enumerate(y):
        ax.annotate(txt, (x[i] - 20, txt + annotate_move_y))
    plt.grid(True)
    ax.set_ylim([0, ax.get_ylim()[1] + annotate_move_y * 2])
    ax.set_ylabel("Odległość Levensteina")
    ax.set_xlabel("Ilość błędów negatywnych[%]")
    plt.savefig("plots/4.png")


def positive_negative_plot():
    """Generates plot by dna size"""
    x = [2, 5, 10, 20]
    y = []
    for l in x:
        dist = 0
        for _ in range(3):
            _, best_solution, original_seq = regenerate_sequence_with_ants(7, 200, 10, 10, negative=l, positive=l)
            dist += td.levenshtein.distance(best_solution[1], original_seq)
        y.append(dist/3)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(x, y)
    annotate_move_y = math.floor(max(y)/20)
    for i, txt in enumerate(y):
        ax.annotate(txt, (x[i] - 20, txt + annotate_move_y))
    plt.grid(True)
    ax.set_ylim([0, ax.get_ylim()[1] + annotate_move_y * 2])
    ax.set_ylabel("Odległość Levensteina")
    ax.set_xlabel("Ilość błędów pozytywnych i negatywnych[%]")
    plt.savefig("plots/5.png")


def evaporate_plot():
    """Generates plot by dna size"""
    x = list(map(lambda x: 0.1 * x * 100, [1, 2, 4, 8]))
    y = []
    for l in [1, 2, 4, 8]:
        dist = 0
        for _ in range(3):
            _, best_solution, original_seq = regenerate_sequence_with_ants(7, 200, 10, 10, evaporate_coeff=l)
            dist += td.levenshtein.distance(best_solution[1], original_seq)
        y.append(dist/3)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(x, y)
    annotate_move_y = math.floor(max(y)/20)
    for i, txt in enumerate(y):
        ax.annotate(txt, (x[i] - 20, txt + annotate_move_y))
    plt.grid(True)
    ax.set_ylim([0, ax.get_ylim()[1] + annotate_move_y * 2])
    ax.set_ylabel("Odległość Levensteina")
    ax.set_xlabel("Ilość odparowanych feromonów co iterację[%]")
    plt.savefig("plots/6.png")


def oligos_plot():
    """Generates plot by oligo length"""
    x = [7, 8, 9, 10]
    y = []
    for l in x:
        dist = 0
        for _ in range(3):
            _, best_solution, original_seq = regenerate_sequence_with_ants(l, 500, 50, 10)
            dist += td.levenshtein.distance(best_solution[1], original_seq)
        y.append(dist/3)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(x, y)
    annotate_move_y = math.floor(max(y)/20)
    for i, txt in enumerate(y):
        ax.annotate(txt, (x[i] - 20, txt + annotate_move_y))
    plt.grid(True)
    ax.set_ylim([0, ax.get_ylim()[1] + annotate_move_y * 2])
    ax.set_ylabel("Odległość Levensteina")
    ax.set_xlabel("Długość oligonukleotydów")
    plt.savefig("plots/7.png")
