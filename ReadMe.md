# relativistic solver hllc with illum part

-base 3d hllc 
-relativistic solver hllc
-illum part
-joint task


#WARNING!!!!
есть ошибка в расчёта излучения для сферы с внутренней границей. 
Часть границы и возможно внутренних ячеек не обрабатывается! (в решении проскакивает внешняя и внутренняя грань)

при слиянии проектов BuildGraph, MakeStruct, Solve (из Make и Solve) были отброшены структуры 
try_trace, связанные с трасировкой сквозб внутреннюю границу (в BuildGraph оставлены).
 Возможно проблема именно здесь. Если при построении graph в лог выводится сообщение try_...,
внимательно проверять Solve!. Следует запустить тест с инициализаций массива коэффициентов интерполяции и без,
если выходные файлы разные, значит расчет излучения неверный!. 
\note для сплошной сетки проьлема не выявлена.

Нет. эффект проявлется и в односвязной области. Временная затычка --- инициализация структуры 
inter_coef по каждому направлению