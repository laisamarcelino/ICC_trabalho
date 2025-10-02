        CC = gcc

    CFLAGS = -O0
    LFLAGS = -lm

      PROG = cgSolver
      MODULES = utils \
	        pcgc \
                $(PROG)
      OBJS = $(addsuffix .o,$(MODULES))
      SRCS = $(addsuffix .c,$(MODULES)) $(addsuffix .h,$(MODULES))

# Lista de arquivos para distribuição
DISTFILES = *.c *.h Makefile LEIAME
DISTDIR = login1-login2

.PHONY: clean purge dist all

%.o: %.c %.h utils.h
	$(CC) -c $(CFLAGS) $<

$(PROG):  $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)

debug:   CFLAGS+=-D__DEBUG__
debug: $(PROG)

clean:
	@echo "Limpando sujeira ....."
	@rm -rf core *~ *.bak

purge: clean
	@echo "Fazendo a faxina ....."
	@rm -f a.out *.o $(PROG)


dist: purge
	@echo "Gerando arquivo de distribuição ($(DISTDIR).tgz) ..."
	@ln -s . $(DISTDIR)
	@tar -chzvf $(DISTDIR).tgz $(addprefix ./$(DISTDIR)/, $(DISTFILES))
	@rm -f $(DISTDIR)
