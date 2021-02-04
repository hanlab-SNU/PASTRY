#!/bin/bash                                                                                             

awk '{
      numline = numline + 1
     } END{
           printf("%8.1f\n", numline);
          }' "$1"

