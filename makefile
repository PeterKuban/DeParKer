CC		= mpicc
SRC_DIR		= src
EXEC_DIR	= bin
SRC		= main.c 
EXEC		= deparker
DUMP_SRC	= dump_counted.c
DUMP_EXEC	= dump_counted
CFLAGS		= -o3 -fopenmp -lm

all: deparker dump_counted

deparker:
	mkdir -p bin
	$(CC) $(CFLAGS) $(SRC_DIR)/$(SRC) $(SRC_DIR)/sort.c $(SRC_DIR)/kmers.c $(SRC_DIR)/partition.c $(SRC_DIR)/defs.h -o $(EXEC_DIR)/$(EXEC)
	
dump_counted:
	$(CC) $(CFLAGS) $(SRC_DIR)/$(DUMP_SRC) -o $(EXEC_DIR)/$(DUMP_EXEC)
	
	
clean:
	rm -f $(EXEC_DIR)/$(EXEC) $(EXEC_DIR)/$(DUMP_EXEC)