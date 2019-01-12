gen_desk(_,0,[]):-!.
gen_desk(0,J,D):- J>0, J1 is J - 1, gen_desk(8, J1, D),!.
gen_desk(I,J,[c(I,J,n) | D]) :- I>0, I1 is I-1, gen_desk(I1, J, D). 


direction([],1):-!.
direction([1],3):-!.
direction([I|_],I1):- I1 is (I+3) mod 8.

isEnd([]):-!.
isEnd([c(_,_,y)|Desk]):-isEnd(Desk).

nCoord(X,Y,X1,Y1, 0):-X<8,Y>1,X1 is X+1, Y1 is Y - 1,!.
nCoord(X,Y,X1,Y,1):-X<8,X1 is X+1,!.
nCoord(X,Y,X1,Y1,2):-X<8,Y<8,X1 is X+1, Y1 is Y+1,!.
nCoord(X,Y,X,Y1,3):-Y<8,Y1 is Y+1,!.
nCoord(X,Y,X1, Y1,4):-X>1,Y<8,X1 is X-1, Y1 is Y+1,!.
nCoord(X,Y,X1, Y,5):-X>1,X1 is X-1,!.
nCoord(X,Y,X1, Y1, 6):-X>1,Y>1,X1 is X-1, Y1 is Y-1,!.
nCoord(X,Y,X,Y1, 7):-Y>1,Y1 is Y-1,!.


broken(1,1,[],_,Desk):-isEnd(Desk).
broken(X,Y,[X1,Y1|Broken],Dir,Desk):-
    direction(Dir,NDir),
    line(X,Y,X1,Y1,NDir,Desk,Desk1,L),
    L>2,
    broken(X1,Y1,Broken,[NDir|Dir],Desk1).

line(X,Y,X0,Y0,Dir,Desk,Desk0,IJ):-
    nCoord(X,Y,X1,Y1,Dir),
    mark(X1,Y1,Desk,Desk1,J),
    line(X1,Y1,X0,Y0,Dir,Desk1,Desk0,I),
    IJ is I + J.

line(X,Y,X,Y,_,Desk,Desk,0).


mark(X,Y,[c(X,Y,n)|Desk],[c(X,Y,y)|Desk],1):-!.
mark(X,Y,[c(X,Y,y)|Desk],[c(X,Y,y)|Desk],0):-!.
mark(X,Y,[Cell|Desk],[Cell|Desk1],I):-mark(X,Y,Desk,Desk1,I).

run :-
       gen_desk(8,8,Desk),
       broken(1,1,Broken,[],Desk),
       format('Путь: ~w~n', [[1,1|Broken]]),
       fail;
       write('Готово.').