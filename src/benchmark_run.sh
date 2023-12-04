#!/usr/bin/env sh

N=(10 50 100 175 250)

for n in "${N[@]}"; do
  echo "--- Running simulations with $n particles"
  sed -i "s@int const N = .*;@int const N = $n;@g" ../src/problems/argon.cpp ../src/problems/argon_cell.cpp
  make

  echo "--- ContainerVector"
  ./argon > "argon_$n.log"

  echo "--- ContainerCells"
  ./argon_cell > "argon_cell_$n.log"
done
