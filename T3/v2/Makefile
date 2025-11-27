# Compilador e flags
CC = gcc
CFLAGS = -O3 -march=native -mavx -fopt-info-vec -DLIKWID_PERFMON -I$(LIKWID_INCLUDE)
LFLAGS = -lm -L$(LIKWID_LIB) -llikwid # -llikwid

# Nome do executável e arquivos do projeto
PROG = cgSolver
MODULES = utils pcgc sislin helpers $(PROG)
OBJS = $(addsuffix .o, $(MODULES))
SRCS = $(addsuffix .c, $(MODULES)) $(addsuffix .h, $(MODULES))

# Arquivos para empacotamento
DISTFILES = *.c *.h Makefile LEIAME
DISTDIR = lmsr21-rrk24

.PHONY: all clean purge dist debug

# all: compila e gera ./cgSolver
all: $(PROG)

# Regra para criar o executável
$(PROG): $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)

# Compilar cada .c individualmente
%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

# Debug mode
debug: CFLAGS += -D__DEBUG__
debug: $(PROG)

# Limpeza dos arquivos gerados
clean:
	@echo "Limpando arquivos temporários..."
	@rm -f a.out *.o *.bak *~ core $(PROG)

purge: clean
	@echo "Removendo executável e temporários..."

# Empacotamento para entrega
dist: purge
	@echo "Gerando arquivo de distribuição ($(DISTDIR).tgz) ..."
	@ln -s . $(DISTDIR)
	@tar -chzvf $(DISTDIR).tgz $(addprefix ./$(DISTDIR)/, $(DISTFILES))
	@rm -f $(DISTDIR)
