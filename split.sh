#!/bin/bash
grep -f $1 -A1 everything | grep -f $2 -A1 | sed '/^--/d' | sed '/^@/d'
