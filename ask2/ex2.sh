#!/bin/bash

# Εργαστήριο εισαγωγής στον παράλληλο υπολογισμό
# Άσκηση 2 / 2021-22
# Ονοματεπώνυμο: Ζωβοϊλης Δημήτριος-Μάριος
# ΑΜ: 19390064
# Το script αυτό κάνει compile και τρέχει
# το αρχείο ex2.c. Το τρέξιμο γίνεται με την χρήση
# heredoc για να μην χρειάζεται input στα scanf

# compile the source code
mpicc -o ex2 ex2.c

# execute our program with a 4x4 mat and 7 procs
mpirun -np 7 ./ex2 << EOF
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

# execute our program with a 5x5 mat and 3 procs
# mpirun -np 3 ./ex2 << EOF
# 5
# 5
# 1
# 2
# 0
# 0
# 1
# 7
# 2
# 3
# 0
# 9
# 0
# 13
# 2
# 0
# 0
# 0
# 0
# 5
# 0
# 1
# 2
# 3
# 4
# 11
# EOF
