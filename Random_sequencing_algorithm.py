from Original_sequence import EMPTY_ELEMENT, Original_sequence
import random

def Random_sequencing_algorithm (start_row,original_sequence: Original_sequence,Recreating_DNA : str,Visited_verticle : list):
    
    while(len(Recreating_DNA)!=len(original_sequence.sequence)):
        
        for column_number,column in enumerate(start_row):

            if len(Recreating_DNA) == len(original_sequence.sequence):
                break

            if (column != EMPTY_ELEMENT and column != 0 and original_sequence.oligo_list[column_number] not in Visited_verticle):
                
                match column:
                    
                    case 1: #jesli jest cos o wadze 1 gdzie jeszcze nie bylam to
                        Visited_verticle, Recreating_DNA = Recreate_sequence(1,Visited_verticle,original_sequence.oligo_list,column_number,Recreating_DNA)#ustawiam nowe parametry po wykonaniu tej funkcji
                        Random_sequencing_algorithm(original_sequence.oligo_matrix[column_number],original_sequence,Recreating_DNA,Visited_verticle) #wywolujemy funkcje ponownie, szukajac sasiadow tego wierzcholka
                                #!!!!!!Ten pierwszy argument to start czyli caly rzad z pustymi polami, 1 2 i 3. Jest unikatowy z zalozenia  ??????
    
                    case 2:
                        Visited_verticle, Recreating_DNA = Recreate_sequence(2,Visited_verticle,original_sequence.oligo_list,column_number,Recreating_DNA)
                        Random_sequencing_algorithm(original_sequence.oligo_matrix[column_number],original_sequence,Recreating_DNA,Visited_verticle)

                    case 3:
                        Visited_verticle, Recreating_DNA = Recreate_sequence(3,Visited_verticle,original_sequence.oligo_list,column_number,Recreating_DNA)
                        Random_sequencing_algorithm(original_sequence.oligo_matrix[column_number],original_sequence,Recreating_DNA,Visited_verticle)

            else:
             pass
                #TO DO MIEJSCE NA PRZECHOWYWANIE ROZWIAZAN Z ALGORYTMU POMOCNICZEGO
                #TO DO WYBRANIE NAJLEPSZEGO ROZWIAZANIA
                #TO DO ZMUSZENIE PROGRAMU ZEBY POSZEDL WYBRANYM ROZWIAZANIEM
                #TO DO CO JESLI WIERZCHOLEK NIE PROWADZI NIGDZIE BO JEST DEAD END (cofamy sie az do momentu w ktorym gdzies mozemy wejsc w co innego
                #lub losuje cokolwiek gdzie sie jeszcze nie bylo 
                # Searching_for_not_visited(0,...)

        
    # if (len(original_sequence.oligo_list)==len(Visited_verticle)):
    #     return Recreating_DNA
    return Recreating_DNA
        

def Recreate_sequence (how_many_nucleotid : int ,Visited_verticle : list,oligo_list,column_number,Recreating_DNA):
    Visited_verticle.append(oligo_list[column_number]) # dodajemy te sekwencje do miejsc w ktorych bylam
    Nucleotid_to_add = oligo_list[column_number]
    Recreating_DNA += (Nucleotid_to_add[-how_many_nucleotid:]) #Dodaje ostatnia zasade azotowa do rozwiazania
    return Visited_verticle, Recreating_DNA



def Searching_for_not_visited (iteration,start_row, original_sequence : Original_sequence, Visited_verticle,Visited_row_nr,lead_to_unvisited,Road,unvisited_option):
      
    #ogarnij czy to nie powinno byc robione BFS i fifo
    unvisited_column = []
    unvisited_option = {} #to ma byc po to aby ta funkcja chodzila inaczej a nie zeby wywolujac ja 10 razy otrzymac 10 takich samych wynikow
    lead_to_unvisited = 0
    Visited_row_nr = []
    current_row=0
    any_connection=[]
    #Connected= [[]]
    Road = [] #do sledzenia trasy np. 1,3,6 oznacza ze z rzedu 1 idziemy do kolumny 3 wchodzimy do rzedu 3, do kolumny 6, wchodzimy do rzedu 6 
    

    for column_nr,column in enumerate(start_row): 
        
        any_connection.clear()
        unvisited_column.clear()
        if (iteration == 10): #zeby nie krazyc w nieskonczonosc tylko np sprawdzic 10 najblizszych
            return lead_to_unvisited,Road,unvisited_option
       
        if (column == 0):
            current_row=column_nr
        if(column != EMPTY_ELEMENT and column != 0):
            any_connection.append(column_nr)
            if (original_sequence.oligo_list[column_nr] not in Visited_verticle):
                lead_to_unvisited+=1
                unvisited_column.append(column_nr) #tutaj zbieram wszystko w czym mnie nie bylo i co potem chce odwiedzic 
    
    unvisited_option[current_row] = unvisited_column #dodaje do slownika
    # Connected.append (current_row,any_connection)   

    Visited_row_nr.append(current_row) #sprawdzilam juz ile w tym rzedzie jest nieodwiedzonych miejsc wiec odhaczam
                   
    if len(unvisited_column) == 0: # w sytuacji w ktorej w tym rzedzie odwiedzone jest wszystko
        return lead_to_unvisited,Road,unvisited_option

        # ZAPYTAC CZY TO NIE ZA DUZO !!!!!!!!!!!!!!!!!!!!!!!
        # if (any_connection not in Visited_row_nr): #szukam czy jest miejsce ktorego jeszcze nie odwiedzilam w ramach tej funkcji
        #     while (go_to in Visited_row_nr): #jesli jest to losuje tak aby nie krazyc w kolko w tej funkcji tylko gdzies wyjsc
        #         go_to = random.choice(any_connection)
        #         Road.append(go_to)
        #         Searching_for_not_visited(go_to,original_sequence,Visited_verticle,Visited_row_nr,lead_to_unvisited,Road,unvisited_option) 
        
        # elif (any_connection in Visited_row_nr): #jesli juz odwiedzilam w ramach tej funkcji wszystkie z any connection
        #     max_length = 0
        #     max_length_index = -1

        #     for x in any_connection: #dla mozliwych opcji
        #         length = len(unvisited_option[x]) #sprawdzam ktora ma najwiecej nieodwiedzonych
        #         if (length>max_length):
        #             max_length=length
        #             max_length_index=x

        #     if (max_length_index == -1): #tutaj to juz ide sobie gdziekolwiek gdzie jest cos nieodwiedzonego, nawet jak sie nie laczy
        #         #z miejscem w ktorym jestem teraz
                
        #         connected_to_sth = -1
        #         for x in unvisited_option:
        #             if (len(unvisited_option[x]) > 0):
        #                 connected_to_sth=x
                        
        #         if (connected_to_sth!=-1):
        #             return lead_to_unvisited,Road,unvisited_option ##WYJSCIE Z FUNKCJI
                
            
        #     elif (max_length_index != -1):
        #         go_to = unvisited_option[max_length_index]
        #         Road.append(go_to) 
        #         Searching_for_not_visited(go_to,original_sequence,Visited_verticle,Visited_row_nr,lead_to_unvisited,Road,unvisited_option) 

    elif (len(unvisited_column)!= 0):
        go_to = unvisited_option[current_row].pop()  #wyciągam z nieodwiedzonego tam gdzie zaraz pójdę
        Road.append(go_to)
        Searching_for_not_visited(iteration + 1,go_to,original_sequence,Visited_verticle,Visited_row_nr,lead_to_unvisited,Road,unvisited_option)
    
    
    return lead_to_unvisited,Road,unvisited_option #unvisited option zeby wywolujac te funkcje iles razy nie isc ciagle tak samo
    