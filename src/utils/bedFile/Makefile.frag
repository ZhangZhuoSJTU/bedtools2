BUILT_OBJECTS += obj/bedFile.o obj/bedFileHits.o

obj/bedFile.o: src/utils/bedFile/bedFile.cpp obj/bedFile.d
	$(CXX_COMPILE)

obj/bedFileHits.o: src/utils/bedFile/bedFileHits.cpp obj/bedFileHits.d
	$(CXX_COMPILE)
