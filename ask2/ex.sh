#!/bin/bash

# Εργαστήριο εισαγωγής στον παράλληλο υπολογισμό
# Άσκηση 2 / 2021-22
# Ονοματεπώνυμο: Ζωβοϊλης Δημήτριος-Μάριος
# ΑΜ: 19390064
# Το script αυτό κάνει compile και τρέχει
# το αρχείο ex.c. Το τρέξιμο γίνεται με την χρήση
# heredoc για να μην χρειάζεται input στα scanf

# compile the source code
mpicc -o ex ex.c

# execute our program with 4 procs
mpirun -np 4 ./ex << EOF
4
5
1
2
0
1
7
2
3
9
0
13
2
0
0
0
5
EOF