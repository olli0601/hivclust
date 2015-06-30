DEST= $(R_LIBS_USER)/hivclust
EXTERNAL_CODE= $(DEST)/ext
CWD= $$(pwd)

all:	
	echo "current directory is "$(CWD)
	echo "destination directory for external code is "$(EXTERNAL_CODE)
	mkdir -p $(EXTERNAL_CODE)
	cp -R ../misc $(DEST)
	@echo "*********************************************** "
	@echo "********* BUILDING standard-RAxML ************* "
	@echo "*********************************************** "	
	make ExaML-raxml	
	@echo "*********************************************** "
	@echo "********* BUILDING ExaML-parsimonator ********* "
	@echo "*********************************************** "	
	make EXAML-parsimonator
	@echo "*********************************************** "
	@echo "********* BUILDING ExaML ********************** "
	@echo "*********************************************** "	
	make examl	
	@echo "*********************************************** "
	@echo "********* BUILDING BEAST/2.1.0 **************** "
	@echo "*********************************************** "	
	make beast2.jar		
	@echo "*********************************************** "
	@echo "********* BUILDING BEAST/2.sampled ancestors ** "
	@echo "*********************************************** "	
	make beast2-SA.jar			
	@echo "*********************************************** "
	@echo "********* BUILDING HIVCLUST R PACKAGE ********* "
	@echo "*********************************************** "			
	MAKEFLAGS="SHLIB_CXXLD=g++ CXX=g++" R CMD SHLIB -o hivc.so *.cpp
	
3seq:
	echo $(CWD)	
	unzip ../inst/3seq_v1.10812_or.zip
	make -C 3seq_v1.10812_or PTFLAGS=-DPVT='\"'$(DEST)'/\"'
	mv $(CWD)/3seq_v1.10812_or/3seq $(EXTERNAL_CODE)/3seq
	rm -rf 3seq_v1.10812_or
	
ExaML-raxml:
	echo $(CWD)	
	git clone https://github.com/stamatak/standard-RAxML.git
	make -C standard-RAxML -f Makefile.SSE3.gcc
	mv $(CWD)/standard-RAxML/raxmlHPC-SSE3 $(EXTERNAL_CODE)/ExaML-raxml
	rm -rf standard-RAxML	
	
EXAML-parsimonator:	
	echo $(CWD)
	git clone https://github.com/stamatak/Parsimonator-1.0.2.git
	make -C Parsimonator-1.0.2 -f Makefile.SSE3.gcc
	mv $(CWD)/Parsimonator-1.0.2/parsimonator-SSE3 $(EXTERNAL_CODE)/ExaML-parsimonator 	
	rm -rf Parsimonator-1.0.2
	
examl:	
	echo $(CWD)	
	tar -xzvf ../inst/ExaML-3.0.14.tar.gz	
	#git clone https://github.com/stamatak/ExaML.git
	make -C ExaML-3.0.14/parser -f Makefile.SSE3.gcc
	mv $(CWD)/ExaML-3.0.14/parser/parse-examl $(EXTERNAL_CODE)/ExaML-parser
	make -C ExaML-3.0.14/examl -f Makefile.SSE3.gcc
	mv $(CWD)/ExaML-3.0.14/examl/examl $(EXTERNAL_CODE)/ExaML-examl
	rm -rf ExaML-3.0.14
	mv $(EXTERNAL_CODE)/ExaML-examl $(EXTERNAL_CODE)/examl
	
beast2.jar:
	echo $(CWD)
	mv 	../inst/beast2.jar $(EXTERNAL_CODE)/beast2.jar
	
beast2-SA.jar:
	echo $(CWD)
	mv 	../inst/beast2-SA.jar $(EXTERNAL_CODE)/beast2-SA.jar	