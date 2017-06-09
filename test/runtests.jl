using NumDiff
using Base.Test

"""Векторозначная неавтономная функция от многих переменных
f(t, x_1, x_2,...,x_n) -> [f1, f2,...,fn]
"""
function f1(t::Float64, x::Vector{Float64})::Vector{Float64}
  return x .^ 2 + t*x
end

"""Частные производные по всем векторным аргументам функции f1"""
function dxf1(t::Float64, x::Vector{Float64})::Matrix{Float64}
  n = length(x)
  return diagm(2.0 .* x) + t * eye(n,n)
end

"""Векторозначная автономная функция от многих переменных
f(x_1, x_2,...,x_n) -> [f1, f2,...,fn]
"""
function f2(x::Vector{Float64})::Vector{Float64}
  return x .^2
end

"""Частные производные по всем векторным аргументам функции f2"""
function dxf2(x::Vector{Float64})::Matrix{Float64}
  return diagm(2.0 .* x)
end

"""Скалярная автономная функция от многих переменных
f(x_1, x_2,...,x_n) -> f::Float64
"""
function f3(x::Vector{Float64})::Float64
  return dot(x, x)
end

"""Частные производные по всем векторным аргументам функции f3"""
function dxf3(x::Vector{Float64})::Vector{Float64}
  return 2.0 .* x
end

"""Скалярная неавтономная функция от многих переменных
f(t, x_1, x_2,...,x_n) -> f::Float64
"""
function f4(t::Float64, x::Vector{Float64})::Float64
  return dot(x, x) + t
end

"""∂f4/∂x1, ∂f4/∂x2,...,∂f4/∂xn"""
function dxf4(t::Float64, x::Vector{Float64})::Vector{Float64}
  return 2.0 .* x
end

"""Векторознаяная функция от одной переменной"""
function f5(x::Float64)::Vector{Float64}
  return [1.0, 2.0, 3.0] * x
end

"""∂f5/∂x1, ∂f5/∂x2,...,∂f5/∂xn"""
function dxf5(x::Float64)::Vector{Float64}
  return [1.0, 2.0, 3.0]
end

# Тестирование
t_0 = 0.5
x_0 = [1.0, 2.0, 3.0]
# для всех порядков
@testset for od in 2:2:8
  # для всех переменных
  @testset for i in 1:length(x_0)
    @test NumDiff.derivative(f1, t_0, x_0; k=i, order=od) ≈ dxf1(t_0, x_0)[:, i]
    @test NumDiff.derivative(f2, x_0; k=i, order=od) ≈ dxf2(x_0)[:, i]
    @test NumDiff.derivative(f3, x_0; k=i, order=od) ≈ dxf3(x_0)[i]
    @test NumDiff.derivative(f4, t_0, x_0; k=i, order=od) ≈ dxf4(t_0, x_0)[i]
  end
  @test NumDiff.derivative(f5, 2.0; order=8) ≈ dxf5(2.0)
end
