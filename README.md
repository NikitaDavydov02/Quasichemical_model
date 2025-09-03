# Quasichemical_model
Программа для численного нахождения профиля распределения полимера в растворе согласно модели Схойтьенса-Флира (10.1021/j100475a012).
Ввод параметров осуществляется через файл input.txt в следующем формате:
20 //number of layers in a cell
1000 //number of segments in a polymer chain
0.5 //flory parameter between polymer and solvent
1 //adsorption parameter
0.01 //bulk polymer volume fraction
0 //output intermediate info to console
1000 //maximum number of iterations
1e-6 //required error
0.01 //newthon break coeffitient

Вывод результатов (объемной доли полимера в каждом из слоев) осуществляеся в консоль, а также в файл output.txt

