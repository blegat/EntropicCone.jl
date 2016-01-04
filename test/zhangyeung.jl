G = polymatroidcone(4)
@test matusrentropy(1, 14) in G
@test matusrentropy(1, 23) in G
@test matusrentropy(1, 24) in G
@test matusrentropy(2, 3) in G
@test matusrentropy(2, 4) in G
@test invalidfentropy(12) in G
