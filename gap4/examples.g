A2:= rec(
  gens:= [1,2],
  invr:= [1,2],
  rels:= [
    [[1,2,1], [2,1,2]],
  ],
  sbgp:= [[2]],
);

A3:= rec(
  gens:= [1,2,3],
  invr:= [1,2,3],
  rels:= [
    [[1,2,1], [2,1,2]],
    [[1,3],[3,1]],
    [[2,3,2],[3,2,3]],
  ],
  sbgp:= [[2],[3]],
);


G4 := rec(
  gens:= [1,2,3,4],
  invr:= [3,4,1,2],
  rels:= [
    [[1,1,1], []],
    [[2,2,2], []],
    [[1,2,1], [2,1,2]],
  ],
  sbgp:= [[1]],
);


G5 := rec(
  gens:= [1,2,3,4],
  invr:= [3,4,1,2],
  rels:= [
    [[1,1,1], []],
    [[2,2,2], []],
    [[1,2,1,2], [2,1,2,1]],
  ],
  sbgp:= [[1]],
);

G6 := rec(
  gens:= [1,2,3],
  invr:= [1,3,2],
  rels:= [
    [[2,2,2], []],
    [[1,2,1,2,1,2], [2,1,2,1,2,1]],
  ],
  sbgp:= [[2]],
);


G7 := rec(
  gens:= [1,2,3,4,5],
  invr:= [1,4,5,2,3],
  rels:= [
    [[2,2,2], []],
    [[3,3,3], []],
    [[1,2,3], [2,3,1]],
    [[2,3,1], [3,1,2]],
  ],
  sbgp:= [[3]],
);


G8 := rec(
  gens:= [1,2,3,4],
  invr:= [3,4,1,2],
  rels:= [
    [[1,1,1,1], []],
    [[2,2,2,2], []],
    [[1,2,1], [2,1,2]],
  ],
  sbgp:= [[1]],
);


G9 := rec(
  gens:= [1,2,3],
  invr:= [1,3,2],
  rels:= [
    [[2,2,2,2], []],
    [[1,2,1,2,1,2], [2,1,2,1,2,1]],
  ],
  sbgp:= [[2]],
);


G10 := rec(
  gens:= [1,2,3,4],
  invr:= [3,4,1,2],
  rels:= [
    [[1,1,1], []],
    [[2,2,2,2], []],
    [[1,2,1,2], [2,1,2,1]],
  ],
  sbgp:= [[2]],
);


G11 := rec(
  gens:= [1,2,3,4,5],
  invr:= [1,4,5,2,3],
  rels:= [
    [[2,2,2], []],
    [[3,3,3,3], []],
    [[1,2,3], [2,3,1]],
    [[2,3,1], [3,1,2]],
  ],
  sbgp:= [[3]],
);


G12 := rec(
  gens:= [1,2,3],
  invr:= [1,2,3],
  rels:= [
    [[1,2,3,1], [2,3,1,2]],
    [[2,3,1,2], [3,1,2,3]],
  ],
  sbgp:= [[1]],
);


G13 := rec(
  gens:= [1,2,3],
  invr:= [1,2,3],
  rels:= [
    [[1,2,3,1,2], [3,1,2,3,1]],
    [[2,3,1,2], [3,1,2,3]],
  ],
  sbgp:= [[1]],
);


G24:= rec(
  gens:= [1,2,3],
  invr:= [1,2,3],
  rels:= [
    [[1,2,1], [2,1,2]],
    [[1,3,1], [3,1,3]],
    [[2,3,2,3], [3,2,3,2]],
#     [[2,3,1,2,3,1,2,3,1], [3,2,3,1,2,3,1,2,3]],
    [[1,2,3,1,2,3,1,-3], [2,3,-2,1,2,3,1,2]],
  ],
  sbgp:= [[2], [3]],
);


G333:= rec(
  gens:= [1,2,3,4],
  invr:= [1,2,3,4],
  rels:= [
    [[1,4,1],[4,1,4]],
    [[2,4,2],[4,2,4]],
    [[3,4,3],[4,3,4]],
    [[3,2],[2,1]],
    [[2,1],[1,3]],
  ],
  sbgp:= [[1], [2], [3]],
);


M12:= rec(
  gens:= [ 1,2,3,4 ],
  invr:= [ 2,1,3,4 ],
  rels:= [
    [[ 1,1,1,1,1,1,1,1,1,1,1 ], []],
    [[ 1,3,1,3,1,3 ], []],
    [[ 1,4,1,4,1,4 ], []],
    [[ 3,4,3,4,3,4,3,4,3,4,3,4,3,4,3,4,3,4,3,4 ], []],
    [[ 1,1,3,4,3,4,1,4,3,4,3 ], []],
  ],
  sbgp:=[[1],[2],[3]],
);
