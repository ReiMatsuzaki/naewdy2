include ../local.mk

# ==== Common ====
FF:=${FF}
MODS0=const write_title err_handle timer sys convert istream fjson math timestep 

# -- Directories --
SRC=${NAEWDY_PATH}/src
BUILD=${NAEWDY_PATH}/build
TESTS=${NAEWDY_PATH}/tests

# -- compile --
%.o: %.f
	${FC} ${FF} -c $< -o $@

%.o: %.f90
	${FC} ${FF} -c $< -o $@

# -- Clean --
clean:
	rm -f *.o
	rm -f *.mod

# ==== NAEWDY ====
# -- config --
MODULES=naewdy
GAMESS_PATH=${NAEWDY_PATH}/gamess/gamess_naewdy
GAMESS_VER=00
TARGET_GAMESS=gamess.${GAMESS_VER}.x

# -- object --
GAMESS_SRCS=$(addprefix ${GAMESS_PATH}/source/, $(addsuffix .src,${MODULES}))

# -- gamess --
${GAMESS_PATH}/source/naewdy.src: ${SRC}/naewdy.f90
	${FC} -E -cpp ${SRC}/naewdy.f90 > ${GAMESS_PATH}/source/naewdy.src

gms:
#	echo "analyse macros and copy source to gamess dirs "
	${FC} -E -cpp ${SRC}/naewdy.f90 > ${GAMESS_PATH}/source/naewdy.src
	echo "start gamess compile"
	make -C ${GAMESS_PATH} gamess
