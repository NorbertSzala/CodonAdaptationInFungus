#!/bin/python3

# Decyduje sie na wyliczenie sredniej dla wszystkich bialek rybosomowych reprezentujacych duzą oraz mala podjednostke czyli 40S oraz 60S poniewaz nie znalazlem takiego bialka z tych podjednostek ktore byloby obecne we wszystkich organizmach (for i in {1..50}; do grep -l "60S ribosomal protein L$i" * | wc -l ; done)

#grep "protein=60S ribosomal protein\|protein=40S ribosomal protein" * | grep -v "mitochondrial" | wc -l # it gives us 8735 ribosomal non-mitochondrial proteins



input_path = ""
output_path = ""