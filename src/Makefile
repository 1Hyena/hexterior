NAME    = hexterior
CC      = g++
PROF    = -O2
ALLEGRO = -L/usr/lib -lallegro -lallegro_font -lallegro_primitives -lallegro_color -lallegro_image
C_FLAGS = -std=c++11 -Wall -Wextra -Werror -Wno-error=unused-variable -pedantic-errors $(PROF)
L_FLAGS = -lm -lstdc++ $(ALLEGRO) $(PROF)
OBJ_DIR = obj

SRC_FILES := $(wildcard *.cpp)
O_FILES   := $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))

OUT = ../$(NAME)

all:
	@$(MAKE) dynamic -s

dynamic: $(O_FILES)
	@printf "\033[1;33mMaking \033[37m   ...."
	$(CC) -o $(OUT) $(O_FILES) $(L_FLAGS)
	@printf "\033[1;32m %s DONE!\033[0m\n" $(NAME)

$(OBJ_DIR)/%.o: %.cpp
		@printf "\033[1m\033[31mCompiling \033[37m....\033[34m %-20s\t\033[33m%6s\033[31m lines\033[0m \n" $*.cpp "`wc -l $*.cpp | cut -f1 -d' '`"
		@$(CC) $< $(C_FLAGS) -c -o $@

clean:
	@printf "\033[1;36mCleaning \033[37m ...."
	@rm -f $(O_FILES) $(OUT) *~ *.bak *.orig *.rej
	@printf "\033[1;37m $(NAME) cleaned!\033[0m\n"
