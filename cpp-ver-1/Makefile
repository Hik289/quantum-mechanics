TQ: main.cpp obj/func.o obj/para.o obj/head.o
	g++ -o TQ main.cpp obj/func.o obj/para.o obj/head.o
obj/func.o: func.cpp
	g++ -c -o obj/func.o func.cpp
obj/para.o: para.cpp
	g++ -c -o obj/para.o para.cpp
obj/head.o: head.cpp
	g++ -c -o obj/head.o head.cpp
clean:
	-rm -f obj/*.o
