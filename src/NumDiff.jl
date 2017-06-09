"""Функции для численного дифференцирования с использованием центральных разностей"""
module NumDiff
"""
Нахождение частной производной первого порядка
от векторозначной функции многих переменных

``derivative(func, [t,] x, k; p, h, order) -> df``

Для вычисления производных используется метод центральных разностей.

# Аргументы

* `func::Function`: функция в виде func(t, x) или func(x).
* `t::Float64`: обциональный скалярный аргумент функции (время).
* `x::Vector{Float64}`: векторный аргумент функции `x`.
* `k::Int`: номер переменной, по которой необходимо взять производную.
* `h::Float64`: шаг приращения, который необходимо использовать.
* `order::Int`: порядок точности {2, 4, 6, 8}.

# Возвращаемые значения

`df::Vector{Float64}`: вектор частных производных.
"""
function derivative(func::Function, t::Float64, x::Vector{Float64}; k=1, h=0.1, order=4)
  # Используются центральные конечные разности коэффициенты этих разностей задаются массивом
  # Заполняем массив коэффициентов
  a = [[-0.5, 0.0, 0.5], # Порядка 2
       [1.0/12.0, -2.0/3.0, 0.0, 2.0/3.0, -1.0/12.0], # Порядка 4
       [-1.0/60.0, 3.0/20.0, -3.0/4.0, 0.0, 3.0/4.0, -3.0/20.0, 1.0/60.0], # Порядка 6
       [1.0/280.0, -4.0/105.0, 1.0/5.0, -4.0/5.0, 0.0, 4.0/5.0, -1.0/5.0, 4.0/105.0, -1.0/280.0]] # Порядка 8

  # Проверяем по возвращаемому значению,
  # является ли функция векторной или скалярной
  if isa(func(t, x), Vector{Float64})
    df = zeros(func(t, x))
  elseif isa(func(t, x), Float64)
    df = 0.0
  end

  Ih = h * eye(Float64, length(x))

  if order == 2
    for j = 0:order
      df = df + a[1][j+1] * func(t, x + (j - 1) * Ih[:, k])
    end
  elseif order == 4
    for j = 0:order
      df = df + a[2][j+1] * func(t, x + (j - 2) * Ih[:, k])
    end
  elseif order == 6
    for j = 0:order
      df = df + a[3][j+1] * func(t, x + (j - 3) * Ih[:, k])
    end
  elseif order == 8
    for j = 0:order
      df = df + a[4][j+1] * func(t, x + (j - 4) * Ih[:, k])
    end
  else
    println(STDERR, "Variable order = {2, 4, 6, 8}! Using order=4")
    for j = 0:order
      df = df + a[2][j+1] * func(t, x + (j - 2) * Ih[:, k])
    end
  end

  return df / h
end

# Частные производные от функции без параметра t
function derivative(func::Function, x::Vector{Float64}; k=1, h=0.1, order=4)
  # Заполняем массив коэффициентов
  a = [[-0.5, 0.0, 0.5], # Порядка 2
       [1.0/12.0, -2.0/3.0, 0.0, 2.0/3.0, -1.0/12.0], # Порядка 4
       [-1.0/60.0, 3.0/20.0, -3.0/4.0, 0.0, 3.0/4.0, -3.0/20.0, 1.0/60.0], # Порядка 6
       [1.0/280.0, -4.0/105.0, 1.0/5.0, -4.0/5.0, 0.0, 4.0/5.0, -1.0/5.0, 4.0/105.0, -1.0/280.0]] # Порядка 8

  # Проверяем по возвращаемому значению,
  # является ли функция векторной или скалярной
  if isa(func(x), Vector{Float64})
    df = zeros(func(x))
  elseif isa(func(x), Float64)
    df = 0.0
  end

  Ih = h * eye(Float64, length(x))

  if order == 2
    for j = 0:order
      df = df + a[1][j+1] * func(x + (j - 1) * Ih[:, k])
    end
  elseif order == 4
    for j = 0:order
      df = df + a[2][j+1] * func(x + (j - 2) * Ih[:, k])
    end
  elseif order == 6
    for j = 0:order
      df = df + a[3][j+1] * func(x + (j - 3) * Ih[:, k])
    end
  elseif order == 8
    for j = 0:order
      df = df + a[4][j+1] * func(x + (j - 4) * Ih[:, k])
    end
  else
    println(STDERR, "Variable order = {2, 4, 6, 8}! Using order=4")
    for j = 0:order
      df = df + a[2][j+1] * func(x + (j - 2) * Ih[:, k])
    end
  end

  return df / h
end

# Частные производные от функции одного аргумента без параметра t
function derivative(func::Function, x::Float64; h=0.1, order=4)
  # Заполняем массив коэффициентов
  a = [[-0.5, 0.0, 0.5], # Порядка 2
       [1.0/12.0, -2.0/3.0, 0.0, 2.0/3.0, -1.0/12.0], # Порядка 4
       [-1.0/60.0, 3.0/20.0, -3.0/4.0, 0.0, 3.0/4.0, -3.0/20.0, 1.0/60.0], # Порядка 6
       [1.0/280.0, -4.0/105.0, 1.0/5.0, -4.0/5.0, 0.0, 4.0/5.0, -1.0/5.0, 4.0/105.0, -1.0/280.0]] # Порядка 8

  # Проверяем по возвращаемому значению,
  # является ли функция векторной или скалярной
  if isa(func(x), Vector{Float64})
    df = zeros(func(x))
  elseif isa(func(x), Float64)
    df = 0.0
  end

  if order == 2
    for j = 0:order
      df = df + a[1][j+1] * func(x + (j - 1) * h)
    end
  elseif order == 4
    for j = 0:order
      df = df + a[2][j+1] * func(x + (j - 2) * h)
    end
  elseif order == 6
    for j = 0:order
      df = df + a[3][j+1] * func(x + (j - 3) * h)
    end
  elseif order == 8
    for j = 0:order
      df = df + a[4][j+1] * func(x + (j - 4) * h)
    end
  else
    println(STDERR, "Variable order = {2, 4, 6, 8}! Using order=4")
    for j = 0:order
      df = df + a[2][j+1] * func(x + (j - 2) * h)
    end
  end

  return df / h
end

"""
Нахождение матрицы Якоби от векторозначной функции многих переменных

``jacobian(func, t, x; h=0.1, order=4) -> J``

# Аргументы

* `func::Function`: функция в виде func(t, x, p).
* `t::Float64`: скалярный аргумент функции (время).
* `x::Vector{Float64}`: векторный аргумент функции `x`.
* `p::Tuple`: параметры функции (даже если их нет).
* `h::Float64`: шаг приращения, который необходимо использовать.
* `order::Int`: порядок точности {2, 4, 6, 8}.

# Возвращаемые значения

* `J::Matrix{Float64}(n, n)`: матрица Якоби.
"""
function jacobian(func::Function, t::Float64, x::Vector{Float64}; h=0.1, order=4)
  n = length(x)
  J = Matrix{Float64}(n, n)

  for i in 1:n
    J[:, i] = derivative(func, t, x; k=i, h=h, order=order)
  end
  return J
end
end # module
