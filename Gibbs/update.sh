#!/bin/sh
#$ -cwd

ctags *.cpp *.C *.h
cscope -b

