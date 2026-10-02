BUILT_OBJECTS += obj/bamToBed.o obj/bamToBedHelp.o obj/bamToBedConvert.o

obj/bamToBed.o: src/bamToBed/bamToBed.cpp obj/bamToBed.d
	$(CXX_COMPILE)

obj/bamToBedHelp.o: src/bamToBed/bamToBedHelp.cpp obj/bamToBedHelp.d
	$(CXX_COMPILE)

obj/bamToBedConvert.o: src/bamToBed/bamToBedConvert.cpp obj/bamToBedConvert.d
	$(CXX_COMPILE)
