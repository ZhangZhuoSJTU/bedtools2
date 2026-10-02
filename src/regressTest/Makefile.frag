BUILT_OBJECTS += obj/regressTestMain.o obj/RegressTest.o obj/RegressTestMemory.o

obj/regressTestMain.o: src/regressTest/regressTestMain.cpp obj/regressTestMain.d
	$(CXX_COMPILE)

obj/RegressTest.o: src/regressTest/RegressTest.cpp obj/RegressTest.d
	$(CXX_COMPILE)

obj/RegressTestMemory.o: src/regressTest/RegressTestMemory.cpp obj/RegressTestMemory.d
	$(CXX_COMPILE)
