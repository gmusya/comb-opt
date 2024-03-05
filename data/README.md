`./tsp-random --input $DATA_DIR/input.txt --output $DATA_DIR/output-random.txt`

`./tsp-spanning --input $DATA_DIR/input.txt --output $DATA_DIR/output-spanning.txt`

`./tsp-local --input $DATA_DIR/input.txt --output $DATA_DIR/output-local.txt --initial_route $DATA_DIR/output-spanning.txt`

`./tsp-iteratively --input $DATA_DIR/input.txt --output $DATA_DIR/output-iteratively.txt`

`./lp --input $DATA_DIR/input.txt --output $DATA_DIR/output-milp.txt`
