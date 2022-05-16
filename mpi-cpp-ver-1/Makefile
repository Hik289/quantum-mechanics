TQ: main.cpp obj/func.o obj/para.o obj/head.o
	mpic++ -o TQ main.cpp obj/func.o obj/para.o obj/head.o
obj/func.o: func.cpp
	mpic++ -c -o obj/func.o func.cpp
obj/para.o: para.cpp
	mpic++ -c -o obj/para.o para.cpp
obj/head.o: head.cpp
	mpic++ -c -o obj/head.o head.cpp
clean:
	-rm -f obj/*.o
