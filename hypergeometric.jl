using Nemo, Test
import Nemo: libarb, libflint, arf_struct, arb_struct, mag_struct
import Nemo: expressify, iszero, isfinite, contains_zero
import Nemo: zero!, one!, add!, sub!, mul!, div!, inv!, addmul!,
       coeff, setcoeff!, mullow, gamma, rgamma, solve, overlaps, sub, derivative,
       rising_factorial, nbits

macro debug()
  return false
end

const ARF_PREC_EXACT = typemax(Int)

function precision(R::AcbField)
  return R.prec
end

#### flint/arb types: arb, acb, arb_poly, acb_poly, and acb_mat ################


mutable struct fmpz_poly_mock
  coeffs::Ptr{Int}
  alloc::Int
  len::Int
end

mutable struct fmpz_mock
  data::Int
end

mutable struct nmag
  rad_exp::Int # fmpz
  rad_man::UInt

  function nmag()
    z = new()
    ccall((:mag_init, libarb), Nothing,
          (Ref{nmag},),
          z)
    finalizer(_mag_clear_fxn, z)
    return z
  end
end

function _mag_clear_fxn(a::nmag)
  ccall((:mag_clear, libarb), Nothing, (Ref{nmag},), a)
end

mutable struct narb
  mid_exp::Int # fmpz
  mid_size::UInt # mp_size_t
  mid_d1::UInt # mantissa_struct
  mid_d2::UInt
  rad_exp::Int # fmpz
  rad_man::UInt

  function narb()
    z = new()
    ccall((:arb_init, libarb), Nothing,
          (Ref{narb}, ),
          z)
    finalizer(_arb_clear_fn, z)
    return z
  end
end

function _arb_clear_fn(x::narb)
  ccall((:arb_clear, libarb), Nothing, (Ref{narb}, ), x)
end

mutable struct nacb
  real_mid_exp::Int # fmpz
  real_mid_size::UInt # mp_size_t
  real_mid_d1::Int # mantissa_struct
  real_mid_d2::Int
  real_rad_exp::Int # fmpz
  real_rad_man::UInt
  imag_mid_exp::Int # fmpz
  imag_mid_size::UInt # mp_size_t
  imag_mid_d1::Int # mantissa_struct
  imag_mid_d2::Int
  imag_rad_exp::Int # fmpz
  imag_rad_man::UInt

  function nacb()
    z = new()
    ccall((:acb_init, libarb), Nothing,
        (Ref{nacb}, ),
        z)
    finalizer(_acb_clear_fn, z)
    return z
  end
end

function _acb_clear_fn(x::nacb)
  ccall((:acb_clear, libarb), Nothing, (Ref{nacb}, ), x)
end


mutable struct nacb_poly
  coeffs::Ptr{Nothing}
  length::Int
  alloc::Int

  function nacb_poly()
    z = new()
    ccall((:acb_poly_init, libarb), Nothing, (Ref{nacb_poly}, ), z)
    finalizer(_acb_poly_clear_fn, z)
    return z
  end
end

function _acb_poly_clear_fn(x::nacb_poly)
  ccall((:acb_poly_clear, libarb), Nothing, (Ref{nacb_poly}, ), x)
end


mutable struct narb_poly
  coeffs::Ptr{Nothing}
  length::Int
  alloc::Int

  function narb_poly()
    z = new()
    ccall((:arb_poly_init, libarb), Nothing, (Ref{narb_poly}, ), z)
    finalizer(_arb_poly_clear_fn, z)
    return z
  end
end

function _arb_poly_clear_fn(x::narb_poly)
  ccall((:arb_poly_clear, libarb), Nothing, (Ref{narb_poly}, ), x)
end


mutable struct narb_mat
  entries::Ptr{Nothing}
  r::Int
  c::Int
  rows::Ptr{Nothing}

  function narb_mat(r::Int, c::Int)
    z = new()
    ccall((:acb_mat_init, libarb), Nothing,
          (Ref{narb_mat}, Int, Int),
          z, r, c)
    finalizer(_arb_mat_clear_fn, z)
    return z
  end
end

function _arb_mat_clear_fn(a::narb_mat)
  ccall((:arb_mat_clear, libarb), Nothing, (Ref{narb_mat}, ), a)
end


mutable struct nacb_mat
  entries::Ptr{Nothing}
  r::Int
  c::Int
  rows::Ptr{Nothing}

  function nacb_mat(r::Int, c::Int)
    z = new()
    ccall((:acb_mat_init, libarb), Nothing,
          (Ref{nacb_mat}, Int, Int),
          z, r, c)
    finalizer(_acb_mat_clear_fn, z)
    return z
  end
end

function _acb_mat_clear_fn(a::nacb_mat)
  ccall((:acb_mat_clear, libarb), Nothing, (Ref{nacb_mat}, ), a)
end

#### constructors ##############################################################

function narb(a::Int)
  z = narb()
  ccall((:arb_set_si, libarb), Nothing,
        (Ref{narb}, Int),
        z, a)
  return z
end

function nacb(a::Int)
  z = nacb()
  ccall((:acb_set_si, libarb), Nothing,
        (Ref{nacb}, Int),
        z, a)
  return z
end

function nacb(a::acb)
  z = nacb()
  ccall((:acb_set, libarb), Nothing,
        (Ref{nacb}, Ref{acb}),
        z, a)
  return z
end

function nacb(a::nacb)
  return a
end

function (R::AcbField)(a::nacb)
  z = R()
  ccall((:acb_set_round, libarb), Nothing,
        (Ref{acb}, Ref{nacb}, Int),
        z, a, precision(R))
  return z
end

function nacb(x::narb)
  z = nacb()
  ccall((:acb_set_arb, libarb), Nothing,
        (Ref{nacb}, Ref{narb}),
        z, x)
  return z
end

function nacb(x::narb, y::narb)
  z = nacb()
  ccall((:acb_set_arb_arb, libarb), Nothing,
        (Ref{nacb}, Ref{narb}, Ref{narb}),
        z, x, y)
  return z
end

function nacb(a::fmpz, p::Int)
  z = nacb()
  ccall((:acb_set_round_fmpz, libarb), Nothing,
        (Ref{nacb}, Ref{fmpz}, Int),
        z, a, p)
  return z  
end

function narb_poly(c0::narb)
  z = narb_poly()
  ccall((:arb_poly_set_coeff_arb, libarb), Nothing,
        (Ref{narb_poly}, Int, Ref{narb}),
        z, 0, c0)
  return z
end

function narb_poly(c0::narb, c1::Int)
  z = narb_poly()
  ccall((:arb_poly_set_coeff_si, libarb), Nothing,
        (Ref{narb_poly}, Int, Int),
        z, 1, c1)
  ccall((:arb_poly_set_coeff_arb, libarb), Nothing,
        (Ref{narb_poly}, Int, Ref{narb}),
        z, 0, c0)
  return z
end

function narb_poly(c::Array{narb})
  z = narb_poly()
  for i in 1:length(c)
    setcoeff!(z, i-1, c[i])
  end
  return z
end

function nacb_poly(c0::nacb, c1::nacb)
  z = nacb_poly()
  ccall((:acb_poly_set_coeff_acb, libarb), Nothing,
        (Ref{nacb_poly}, Int, Ref{nacb}),
        z, 1, c1)
  ccall((:acb_poly_set_coeff_acb, libarb), Nothing,
        (Ref{nacb_poly}, Int, Ref{nacb}),
        z, 0, c0)
  return z
end

function nacb_poly(c0::nacb, c1::Int)
  z = nacb_poly()
  ccall((:acb_poly_set_coeff_si, libarb), Nothing,
        (Ref{nacb_poly}, Int, Int),
        z, 1, c1)
  ccall((:acb_poly_set_coeff_acb, libarb), Nothing,
        (Ref{nacb_poly}, Int, Ref{nacb}),
        z, 0, c0)
  return z
end

function nacb_poly(c0::Int)
  z = nacb_poly()
  ccall((:acb_poly_set_coeff_si, libarb), Nothing,
        (Ref{nacb_poly}, Int, Int),
        z, 0, c0)
  return z
end

function nacb_poly(c0::Int, c1::Int)
  z = nacb_poly()
  ccall((:acb_poly_set_coeff_si, libarb), Nothing,
        (Ref{nacb_poly}, Int, Int),
        z, 1, c1)
  ccall((:acb_poly_set_coeff_acb, libarb), Nothing,
        (Ref{nacb_poly}, Int, Ref{nacb}),
        z, 0, c0)
  return z
end

function nacb_poly(c::Array{nacb})
  z = nacb_poly()
  for i in 1:length(c)
    setcoeff!(z, i-1, c[i])
  end
  return z
end

#### output ####################################################################


function set!(z::narb, a::AbstractString, p::Int)
  s = string(a)
  err = ccall((:arb_set_str, libarb), Cint,
              (Ref{narb}, Ptr{UInt8}, Int),
              z, s, p)
  err == 0 || error("Invalid real string: $(repr(s))")
  return z
end

function native_string(x::narb)
  d = precision(x)
  d = clamp(d, 10, 200)
  d = trunc(Int, d*0.30103)
  cstr = ccall((:arb_get_str, Nemo.libarb), Ptr{UInt8},
         (Ref{narb}, Int, UInt),
         x, Int(d), UInt(0))
  r = unsafe_string(cstr)
  ccall((:flint_free, libflint), Nothing, (Ptr{UInt8},), cstr)
  return r
end

function Base.show(io::IO, a::narb)
  print(io, native_string(a))
end

function Base.show(io::IO, a::nacb)
  print(io, native_string(real(a)))
  print(io, " + i*")
  print(io, native_string(imag(a)))
end

function Base.show(io::IO, a::Union{narb_poly, nacb_poly})
  first = true
  for i in 0:length(a)-1
    if !first
      print(io, " + ")
    end
    first = false
    print(io, "(")
    print(io, coeff(a, i))
    print(io, ")*ε^"*string(i))
  end
  if first
    print(io, "0")
  end
end

function Base.show(io::IO, a::Union{narb_mat,nacb_mat})
  println(io, string(nrows(a))*" by "*string(ncols(a)))
  println(io, " [")
  for i in 1:nrows(a), j in 1:ncols(a)
    println(io, a[i,j])
  end
  println(io, "]")
end

Nemo.AbstractAlgebra.needs_parentheses(x::Union{narb, nacb, nacb_poly, nacb_mat}) = true


#### boring arithmetic #########################################################

function radius(a::narb)
  GC.@preserve a begin
    r = ccall((:arb_rad_ptr, libarb), Ptr{mag_struct},
              (Ref{narb}, ),
              a)
    return ccall((:mag_get_d, libarb), Float64,
                 (Ptr{mag_struct}, ),
                 r)
  end
end

function radius(a::nacb)
  return max(radius(real(a)), radius(imag(a)))
end


nrows(a::narb_mat) = a.r
ncols(a::narb_mat) = a.c
nrows(a::nacb_mat) = a.r
ncols(a::nacb_mat) = a.c

function getindex!(z::narb, a::narb_mat, i::Int, j::Int)
  GC.@preserve a begin
    aij = ccall((:arb_mat_entry_ptr, libarb), Ptr{nacb},
                (Ref{narb_mat}, Int, Int),
                a, i - 1, j - 1)
    ccall((:arb_set, libarb), Nothing,
          (Ref{narb}, Ptr{narb}),
          z, aij)
  end
  return z
end

function Base.getindex(a::narb_mat, i::Int, j::Int)
  @assert 0 < i <= nrows(a)
  @assert 0 < j <= ncols(a)
  z = narb()
  getindex!(z, a, i, j)
  return z
end


function getindex!(z::nacb, a::nacb_mat, i::Int, j::Int)
  GC.@preserve a begin
    aij = ccall((:acb_mat_entry_ptr, libarb), Ptr{nacb},
                (Ref{nacb_mat}, Int, Int),
                a, i - 1, j - 1)
    ccall((:acb_set, libarb), Nothing,
          (Ref{nacb}, Ptr{nacb}),
          z, aij)
  end
  return z
end

function Base.getindex(a::nacb_mat, i::Int, j::Int)
  @assert 0 < i <= nrows(a)
  @assert 0 < j <= ncols(a)
  z = nacb()
  getindex!(z, a, i, j)
  return z
end

function Base.setindex!(a::narb_mat, b::narb, i::Int, j::Int)
  GC.@preserve a begin
    aij = ccall((:arb_mat_entry_ptr, libarb), Ptr{narb},
                (Ref{narb_mat}, Int, Int),
                 a, i - 1, j - 1)
    ccall((:arb_set, libarb), Nothing,
          (Ptr{narb}, Ref{narb}),
          aij, b)
  end
  return b
end

function Base.setindex!(a::nacb_mat, b::narb, i::Int, j::Int)
  GC.@preserve a begin
    aij = ccall((:acb_mat_entry_ptr, libarb), Ptr{nacb},
                (Ref{nacb_mat}, Int, Int),
                 a, i - 1, j - 1)
    ccall((:acb_set_arb, libarb), Nothing,
          (Ptr{nacb}, Ref{narb}),
          aij, b)
  end
  return b
end

function Base.setindex!(a::nacb_mat, b::nacb, i::Int, j::Int)
  GC.@preserve a begin
    aij = ccall((:acb_mat_entry_ptr, libarb), Ptr{nacb},
                (Ref{nacb_mat}, Int, Int),
                 a, i - 1, j - 1)
    ccall((:acb_set, libarb), Nothing,
          (Ptr{nacb}, Ref{nacb}),
          aij, b)
  end
  return b
end

function Base.length(a::Union{narb_poly, nacb_poly})
  return a.length
end

function coeff(a::nacb_poly, n::Int)
  @assert n >= 0
  z = nacb()
  ccall((:acb_poly_get_coeff_acb, libarb), Nothing,
        (Ref{nacb}, Ref{nacb_poly}, Int),
        z, a, n)
  return z
end

function coeff(a::narb_poly, n::Int)
  @assert n >= 0
  z = narb()
  ccall((:arb_poly_get_coeff_arb, libarb), Nothing,
        (Ref{narb}, Ref{narb_poly}, Int),
        z, a, n)
  return z
end

function setcoeff!(a::nacb_poly, n::Int, b::nacb)
  @assert n >= 0
  ccall((:acb_poly_set_coeff_acb, libarb), Nothing,
        (Ref{nacb_poly}, Int, Ref{nacb}),
        a, n, b)
  return a
end

function setcoeff!(a::narb_poly, n::Int, b::narb)
  @assert n >= 0
  ccall((:arb_poly_set_coeff_arb, libarb), Nothing,
        (Ref{narb_poly}, Int, Ref{narb}),
        a, n, b)
  return a
end



function precision(x::narb)
  return ccall((:arb_rel_accuracy_bits, libarb), Int,
               (Ref{narb}, ),
               x)
end

function precision(z::nacb)
  return ccall((:acb_rel_accuracy_bits, libarb), Int,
               (Ref{nacb}, ),
               z)
end

function precision_inc(a::Union{narb, nacb}, b::Int)
  p = precision(a)
  return clamp(p, 1, ARF_PREC_EXACT - b) + b
end

function Base.real(a::nacb)
  z = narb()
  ccall((:acb_get_real, libarb), Nothing,
        (Ref{narb}, Ref{nacb}),
        z, a)
  return z
end

function Base.imag(a::nacb)
  z = narb()
  ccall((:acb_get_imag, libarb), Nothing,
        (Ref{narb}, Ref{nacb}),
        z, a)
  return z
end


function convert(::Type{Float64}, x::narb)
  GC.@preserve x begin
    t = ccall((:arb_mid_ptr, libarb), Ptr{arf_struct}, (Ref{narb}, ), x)
    return ccall((:arf_get_d, libarb), Float64, (Ptr{arf_struct}, Int), t, 4)
  end
end

function convert(::Type{ComplexF64}, a::nacb)
  GC.@preserve a begin
    r = ccall((:acb_real_ptr, libarb), Ptr{narb}, (Ref{nacb}, ), a)
    i = ccall((:acb_imag_ptr, libarb), Ptr{narb}, (Ref{nacb}, ), a)
    t = ccall((:arb_mid_ptr, libarb), Ptr{arf_struct}, (Ptr{narb}, ), r)
    u = ccall((:arb_mid_ptr, libarb), Ptr{arf_struct}, (Ptr{narb}, ), i)
    v = ccall((:arf_get_d, libarb), Float64, (Ptr{arf_struct}, Int), t, 4)
    w = ccall((:arf_get_d, libarb), Float64, (Ptr{arf_struct}, Int), u, 4)
  end
  return complex(v, w)
end


function numerator!(z::fmpz, x::fmpq)
  ccall((:fmpq_numerator, libflint), Nothing,
        (Ref{fmpz}, Ref{fmpq}),
        z, x)
  return z
end

function denominator!(z::fmpz, x::fmpq)
  ccall((:fmpq_denominator, libflint), Nothing,
        (Ref{fmpz}, Ref{fmpq}),
        z, x)
  return z
end

function Base.zero(::Type{nacb})
  return nacb()
end

function Base.zero(::Type{narb})
  return narb()
end

function Base.zero(::Type{nacb_poly})
  return nacb_poly()
end

function Base.zero(::Type{narb_poly})
  return narb_poly()
end

function zero!(z::fmpq)
  ccall((:fmpq_set_si, libflint), Nothing,
        (Ref{fmpq}, Int, UInt),
        z, Int(1), UInt(1))
  return z
end

function zero!(z::narb)
  ccall((:arb_zero, libarb), Nothing,
        (Ref{narb}, ),
        z)
  return z
end

function zero!(z::nacb)
  ccall((:acb_zero, libarb), Nothing,
        (Ref{nacb}, ),
        z)
  return z
end

function zero!(z::narb_poly)
  ccall((:arb_zero, libarb), Nothing,
        (Ref{narb_poly}, ),
        z)
  return z
end

function zero!(z::nacb_poly)
  ccall((:acb_poly_zero, libarb), Nothing,
        (Ref{nacb_poly}, ),
        z)
  return z
end

function iszero(a::narb)
  r = ccall((:arb_is_zero, libarb), Cint,
            (Ref{narb}, ),
            a)
  return r != 0
end

function iszero(a::nacb)
  r = ccall((:acb_is_zero, libarb), Cint,
            (Ref{nacb}, ),
            a)
  return r != 0
end

function isnegative(a::narb)
  return 0 != ccall((:arb_is_negative, libarb), Nothing,
                    (Ref{narb},),
                    a)
end

function ispositive(a::narb)
  return 0 != ccall((:arb_is_positive, libarb), Nothing,
                    (Ref{narb},),
                    a)
end

function pos_inf!(z::narb)
  ccall((:arb_pos_inf, libarb), Nothing,
        (Ref{narb},),
        z)
  return z
end

function contains_zero(a::narb)
  r = ccall((:arb_contains_zero, libarb), Cint,
            (Ref{narb}, ),
            a)
  return r != 0
end

function contains_zero(a::nacb)
  r = ccall((:acb_contains_zero, libarb), Cint,
            (Ref{nacb}, ),
            a)
  return r != 0
end

function indeterminate!(z::narb)
  ccall((:arb_indeterminate, libarb), Nothing,
        (Ref{narb}, ),
        z)
  return z
end

function indeterminate!(z::nacb)
  ccall((:acb_indeterminate, libarb), Nothing,
        (Ref{nacb}, ),
        z)
  return z
end

function isfinite(a::nacb)
  r = ccall((:acb_is_finite, libarb), Cint,
            (Ref{nacb}, ),
            a)
  return r != 0
end

function Base.one(::Type{narb})
  return one!(narb())
end

function Base.one(::Type{nacb})
  return one!(nacb())
end

function Base.one(::Type{narb_poly})
  z = narb_poly()
  ccall((:arb_poly_set_si, libarb), Nothing,
        (Ref{narb_poly}, Int),
        z, 1)
  return z  
end

function Base.one(::Type{nacb_poly})
  z = nacb_poly()
  ccall((:acb_poly_set_si, libarb), Nothing,
        (Ref{nacb_poly}, Int),
        z, 1)
  return z  
end

function one!(z::fmpq)
  ccall((:fmpq_set_si, libflint), Nothing,
        (Ref{fmpq}, Int, UInt),
        z, Int(1), UInt(1))
  return z
end

function one!(z::narb)
  ccall((:arb_set_ui, libarb), Nothing,
        (Ref{narb}, UInt),
        z, UInt(1))
  return z
end

function one!(z::nacb)
  ccall((:acb_set_ui, libarb), Nothing,
        (Ref{nacb}, UInt),
        z, UInt(1))
  return z
end


function set!(z::narb, x::fmpq, p::Int)
  ccall((:arb_set_fmpq, libarb), Nothing,
        (Ref{narb}, Ref{fmpq}, Int),
        z, x, p)
  return z
end

function Base.max(x::narb, y::narb, p::Int)
  z = narb()
  ccall((:arb_max, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Ref{narb}, Int),
        z, x, y, p)
  return z
end


function add!(z::fmpq, a::fmpq, b::Int)
  ccall((:fmpq_add_si, libflint), Nothing,
        (Ref{fmpq}, Ref{fmpq}, Int),
        z, a, b)
end

function divexact!(z::fmpz, a::fmpz, b::fmpz)
  @assert divides(a, b)[1]
  ccall((:fmpz_divexact, libflint), Nothing,
        (Ref{fmpz}, Ref{fmpz}, Ref{fmpz}),
        z, a, b)
  return z
end

function mul!(z::fmpz, a::fmpz, b::Int)
  ccall((:fmpz_mul_si, libflint), Nothing,
        (Ref{fmpz}, Ref{fmpz}, Int),
        z, a, b)
  return z
end

function mul!(z::fmpq, a::fmpq, b::Int)
  ccall((:fmpq_mul_si, libflint), Nothing,
        (Ref{fmpq}, Ref{fmpq}, Int),
        z, a, b)
  return z
end

function mul!(z::fmpq, a::fmpq, b::fmpq)
  ccall((:fmpq_mul, libflint), Nothing,
        (Ref{fmpq}, Ref{fmpq}, Ref{fmpq}),
        z, a, b)
  return z
end

function mul!(z::narb, x::narb, y::Int, p::Int)
  ccall((:arb_mul_si, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Int, Int),
        z, x, y, p)
  return z
end

function mul!(z::narb, x::narb, y::fmpz, p::Int)
  ccall((:arb_mul_fmpz, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Ref{fmpz}, Int),
        z, x, y, p)
  return z
end

function mul(a::Union{narb, nacb}, b::Union{fmpz, Int}, p::Int)
  return mul!(typeof(a)(), a, b, p)
end

function mul!(z::nacb, x::nacb, y::fmpz, p::Int)
  ccall((:acb_mul_fmpz, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Ref{fmpz}, Int),
        z, x, y, p)
  return z
end

function mul!(z::nacb, a::nacb, b::narb, p::Int)
  ccall((:acb_mul_arb, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Ref{narb}, Int),
        z, a, b, p)
  return z
end


function mul!(z::narb, x::narb, y::narb, p::Int)
  ccall((:arb_mul, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Ref{narb}, Int),
        z, x, y, p)
  return z
end

function mul!(z::nacb, x::nacb, y::nacb, p::Int)
  ccall((:acb_mul, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Ref{nacb}, Int),
        z, x, y, p)
  return z
end

function mul!(z::nacb_poly, x::nacb_poly, y::nacb, p::Int)
  ccall((:acb_poly_scalar_mul, libarb), Nothing,
        (Ref{nacb_poly}, Ref{nacb_poly}, Ref{nacb}, Int),
        z, x, y, p)
  return z
end

function addmul!(z::narb, a::narb, b::narb, p::Int)
  ccall((:arb_addmul, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Ref{narb}, Int),
        z, a, b, p)
  return z
end

function addmul!(z::nacb, a::nacb, b::nacb, p::Int)
  ccall((:acb_addmul, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Ref{nacb}, Int),
        z, a, b, p)
  return z
end

function submul!(z::narb, a::narb, b::narb, p::Int)
  ccall((:arb_submul, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Ref{narb}, Int),
        z, a, b, p)
  return z
end

function submul!(z::nacb, a::nacb, b::nacb, p::Int)
  ccall((:acb_submul, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Ref{nacb}, Int),
        z, a, b, p)
  return z
end


function div!(z::fmpq, x::fmpq, y::fmpq)
  ccall((:fmpq_div, libflint), Nothing,
        (Ref{fmpq}, Ref{fmpq}, Ref{fmpq}),
        z, x, y)
  return z
end

function div!(z::narb, x::narb, y::Int, p::Int)
  ccall((:arb_div_si, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Int, Int),
        z, x, y, p)
  return z
end

function div!(z::narb, x::narb, y::fmpz, p::Int)
  ccall((:arb_div_fmpz, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Ref{fmpz}, Int),
        z, x, y, p)
  return z
end

function div!(z::narb, a::narb, b::narb, p::Int)
  ccall((:arb_div, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Ref{narb}, Int),
        z, a, b, p)
  return z
end

function div!(z::nacb, a::nacb, b::fmpz, p::Int)
  ccall((:acb_div_fmpz, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Ref{fmpz}, Int),
        z, a, b, p)
  return z
end

function div!(z::nacb, a::nacb, b::narb, p::Int)
  ccall((:acb_div_arb, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Ref{narb}, Int),
        z, a, b, p)
  return z
end

function Base.div(a::T, b::Union{narb, fmpz}, p::Int) where T <: Union{nacb, narb}
  return div!(T(), a, b, p)
end

function mul!(z::T, x::T, y::fmpq, p::Int) where T <: Union{narb, nacb}
  t = fmpz()
  numerator!(t, y)
  mul!(z, x, t, p)
  denominator!(t, y)
  div!(z, z, t, p)
  return z
end

function add!(z::T, x::T, y::fmpq, p::Int) where T <: Union{narb, nacb}
  add!(z, x, T(y, p), p)
end


function mullow(a::narb_poly, b::narb_poly, ord::Int, p::Int)
  z = narb_poly()
  ccall((:arb_poly_mullow, libarb), Nothing,
        (Ref{narb_poly}, Ref{narb_poly}, Ref{narb_poly}, Int, Int),
        z, a, b, ord, p)
  return z  
end

function mullow!(z::nacb_poly, a::nacb_poly, b::nacb_poly, ord::Int, p::Int)
  ccall((:acb_poly_mullow, libarb), Nothing,
        (Ref{nacb_poly}, Ref{nacb_poly}, Ref{nacb_poly}, Int, Int),
        z, a, b, ord, p)
  return z  
end

function mullow(a::nacb_poly, b::nacb_poly, ord::Int, p::Int)
  return mullow!(nacb_poly(), a, b, ord, p)
end

function div_series(a::nacb_poly, b::nacb_poly, ord::Int, p::Int)
  z = nacb_poly()
  ccall((:acb_poly_div_series, libarb), Nothing,
        (Ref{nacb_poly}, Ref{nacb_poly}, Ref{nacb_poly}, Int, Int),
        z, a, b, ord, p)
  return z  
end

function add!(z::nacb_poly, a::nacb_poly, b::nacb_poly, p::Int)
  ccall((:acb_poly_add, libarb), Nothing,
        (Ref{nacb_poly}, Ref{nacb_poly}, Ref{nacb_poly}, Int),
        z, a, b, p)
  return z  
end

function add!(z::narb, x::narb, y::Int, p::Int)
  ccall((:arb_add_si, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Int, Int),
        z, x, y, p)
  return z
end

function add(a::narb, b::Int, p::Int)
  return add!(narb(), a, b, p)
end

function max!(z::narb, x::narb, y::narb, p::Int)
  ccall((:arb_max, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Ref{narb}, Int),
        z, x, y, p)
  return z
end

function add!(z::narb, x::narb, y::narb, p::Int)
  ccall((:arb_add, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Ref{narb}, Int),
        z, x, y, p)
  return z
end

function add(a::narb, b::narb, p::Int)
  return add!(narb(), a, b, p)
end

function sub!(z::narb, x::narb, y::narb, p::Int)
  ccall((:arb_sub, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Ref{narb}, Int),
        z, x, y, p)
  return z
end

function sub(a::narb, b::narb, p::Int)
  return sub!(narb(), a, b, p)
end

function sub(a::nacb, b::nacb, p::Int)
  return sub!(nacb(), a, b, p)
end

function sub(a::Int, b::T, p::Int) where T <: Union{narb, nacb}
  return sub!(T(), a, b, p)
end

function add(a::Int, b::nacb, p::Int)
  return add!(nacb(), b, a, p)
end

function add(b::nacb, a::Int, p::Int)
  return add!(nacb(), b, a, p)
end



function add!(z::narb_poly, x::narb_poly, y::narb_poly, p::Int)
  ccall((:arb_poly_add, libarb), Nothing,
        (Ref{narb_poly}, Ref{narb_poly}, Ref{narb_poly}, Int),
        z, x, y, p)
  return z
end

function add(a::narb_poly, b::narb_poly, p::Int)
  return add!(narb_poly(), a, b, p)
end



function add!(z::nacb, a::nacb, b::Int, p::Int)
  ccall((:acb_add_si, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Int, Int),
        z, a, b, p)
  return z
end

function add!(z::nacb, a::nacb, b::narb, p::Int)
  ccall((:acb_add_arb, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Ref{narb}, Int),
        z, a, b, p)
  return z
end


function add(a::nacb, b::Union{narb, nacb}, p::Int)
  return add!(nacb(), a, b, p)
end

function add!(z::nacb, x::nacb, y::nacb, p::Int)
  ccall((:acb_add, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Ref{nacb}, Int),
        z, x, y, p)
  return z
end

function sub!(z::fmpq, x::fmpq, y::fmpq)
  ccall((:fmpq_sub, libflint), Nothing,
        (Ref{fmpq}, Ref{fmpq}, Ref{fmpq}),
        z, x, y)
  return z
end

function sub!(z::narb, x::narb, y::Int, p::Int)
  ccall((:arb_sub_si, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Int, Int),
        z, x, y, p)
  return z
end

function sub(a::narb, b::Int, p::Int)
  return sub!(narb(), a, b, p)
end

function sub!(z::nacb, x::nacb, y::Int, p::Int)
  ccall((:acb_sub_si, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Int, Int),
        z, x, y, p)
  return z
end

function sub!(z::nacb, a::nacb, b::nacb, p::Int)
  ccall((:acb_sub, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Ref{nacb}, Int),
        z, a, b, p)
  return z
end

function add!(z::narb_poly, a::narb_poly, b::Int, p::Int)
  ccall((:arb_poly_add_si, libarb), Nothing,
        (Ref{narb_poly}, Ref{narb_poly}, Int, Int),
        z, a, b, p)
  return z
end

function add!(z::nacb_poly, a::nacb_poly, b::Int, p::Int)
  ccall((:acb_poly_add_si, libarb), Nothing,
        (Ref{nacb_poly}, Ref{nacb_poly}, Int, Int),
        z, a, b, p)
  return z
end

function add(a::narb_poly, b::Int, p::Int)
  return add!(narb_poly(), a, b, p)
end

function add(a::nacb_poly, b::Union{Int, nacb_poly}, p::Int)
  return add!(nacb_poly(), a, b, p)
end

function add(b::Int, a::narb_poly, p::Int)
  return add!(narb_poly(), a, b, p)
end


function sub!(z::narb_poly, a::Int, b::narb_poly, p::Int)
  return neg!(add!(z, b, Base.checked_neg(a), p)) # TODO
end

function sub(a::Int, b::narb_poly, p::Int)
  return sub!(narb_poly(), a, b, p)
end

function neg!(z::narb, a::narb)
  ccall((:arb_neg, libarb), Nothing,
        (Ref{narb}, Ref{narb}),
        z, a)
  return z
end

function neg!(z::nacb, a::nacb)
  ccall((:acb_neg, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}),
        z, a)
  return z
end

function neg!(z::narb_poly, a::narb_poly)
  ccall((:arb_poly_neg, libarb), Nothing,
        (Ref{narb_poly}, Ref{narb_poly}),
        z, a)
  return z
end

function neg!(z::nacb_poly, a::nacb_poly)
  ccall((:acb_poly_neg, libarb), Nothing,
        (Ref{nacb_poly}, Ref{nacb_poly}),
        z, a)
  return z
end

function neg!(z::Union{narb, nacb, narb_poly, nacb_poly})
  return neg!(z, z)
end

function neg(z::Union{narb, nacb, narb_poly, nacb_poly})
  return neg!(typeof(z)(), z)
end

function sub!(z::nacb, x::Int, y::nacb, p::Int)
  sub!(z, y, x, p)
  neg!(z, z)
  return z
end

function div!(z::nacb, x::nacb, y::nacb, p::Int)
  ccall((:acb_div, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Ref{nacb}, Int),
        z, x, y, p)
  return z
end

function div!(z::narb, x::narb, y::Int, p::Int)
  ccall((:arb_div_si, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Int, Int),
        z, x, y, p)
  return z
end

function div!(z::nacb, x::nacb, y::Int, p::Int)
  ccall((:acb_div_si, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Int, Int),
        z, x, y, p)
  return z
end

function inv!(z::narb, x::narb, p::Int)
  ccall((:arb_inv, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Int),
        z, x, p)
  return z
end

function inv!(z::nacb, x::nacb, p::Int)
  ccall((:acb_inv, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Int),
        z, x, p)
  return z
end

function Base.inv(a::narb, p::Int)
  return inv!(narb(), a, p)
end

function Base.inv(a::nacb, p::Int)
  return inv!(nacb(), a, p)
end

function Base.abs(a::nacb, p::Int)
  return abs!(narb(), a, p)
end

function sqrt!(z::nacb, x::nacb, p::Int)
  ccall((:acb_sqrt, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Int),
        z, x, p)
  return z
end

function abs!(z::narb, a::narb)
  ccall((:arb_abs, libarb), Nothing,
        (Ref{narb}, Ref{narb}),
        z, a)
  return z
end

function abs!(z::narb)
  return abs!(z, z)
end

function abs!(z::narb, a::nacb, p::Int)
  ccall((:acb_abs, libarb), Nothing,
        (Ref{narb}, Ref{nacb}, Int),
        z, a, p)
  return z
end

function log!(z::narb, a::narb, p::Int)
  ccall((:arb_log, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Int),
        z, a, p)
  return z
end

function log!(z::nacb, a::nacb, p::Int)
  ccall((:acb_log, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Int),
        z, a, p)
  return z
end

function pow!(z::narb, a::narb, b::UInt, p::Int)
  ccall((:arb_pow_ui, libarb), Nothing,
      (Ref{narb}, Ref{narb}, UInt, Int),
      z, a, b, p)
  return z
end

function pow!(z::narb, a::narb, b::fmpq, p::Int)
  ccall((:arb_pow_fmpq, libarb), Nothing,
      (Ref{narb}, Ref{narb}, Ref{fmpq}, Int),
      z, a, b, p)
  return z
end

function pow!(z::nacb, a::nacb, b::Int, p::Int)
  ccall((:acb_pow_si, libarb), Nothing,
      (Ref{nacb}, Ref{nacb}, Int, Int),
      z, a, b, p)
  return z
end

function pow!(z::nacb, a::nacb, b::fmpq, p::Int)
  B = narb(b, p)
  ccall((:acb_pow_arb, libarb), Nothing,
      (Ref{nacb}, Ref{nacb}, Ref{narb}, Int),
      z, a, B, p)
  return z
end


function ldexp!(z::nacb, x::nacb, y::Int)
  ccall((:acb_mul_2exp_si, libarb), Nothing,
      (Ref{nacb}, Ref{nacb}, Int),
      z, x, y)
  return z
end

function mul!(z::nacb, a::nacb, b::Int, p::Int)
  ccall((:acb_mul_si, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Int, Int),
        z, a, b, p)
  return z
end

function pow!(z::nacb, a::nacb, b::Int, p::Int)
  ccall((:acb_pow_si, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Int, Int),
        z, a, b, p)
  return z
end

function pow!(z::nacb, a::nacb, b::UInt, p::Int)
  ccall((:acb_pow_ui, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, UInt, Int),
        z, a, b, p)
  return z
end

function mul(x::narb, y::Union{fmpq, narb}, p::Int)
  return mul!(narb(), x, y, p)
end

function mul(x::nacb, y::Union{fmpq, narb, nacb}, p::Int)
  return mul!(nacb(), x, y, p)
end

function Base.sqrt(x::nacb, p::Int)
  return sqrt!(nacb(), x, p)
end

function Base.log(x::nacb, p::Int)
  return log!(nacb(), x, p)
end

function Base.div(x::narb, y::Int, p::Int)
  return div!(narb(), x, y, p)
end

function Base.div(x::nacb, y::nacb, p::Int)
  return div!(nacb(), x, y, p)
end

function pow(x::nacb, y::Union{Int, UInt, fmpz, fmpq, narb, nacb}, p::Int)
  return pow!(nacb(), x, y, p)
end

function pow(a::narb, b::UInt, p::Int)
  return pow!(narb(), a, b, p)
end

function solve(a::nacb_mat, b::nacb_mat, p::Int)
  z = nacb_mat(ncols(a), ncols(b))
  ok = ccall((:acb_mat_solve, libarb), Cint,
             (Ref{nacb_mat}, Ref{nacb_mat}, Ref{nacb_mat}, Int),
             z, a, b, p)
  @assert ok != 0
  return z
end

function add!(z::narb_mat, a::narb_mat, b::narb_mat, p::Int)
  ccall((:arb_mat_add, libarb), Cint,
        (Ref{narb_mat}, Ref{narb_mat}, Ref{narb_mat}, Int),
        z, a, b, p)
  return z
end

function mul!(z::narb_mat, a::narb_mat, b::narb_mat, p::Int)
  ccall((:arb_mat_mul, libarb), Cint,
        (Ref{narb_mat}, Ref{narb_mat}, Ref{narb_mat}, Int),
        z, a, b, p)
  return z
end


function mul(a::nacb_mat, b::nacb_mat, p::Int)
  z = nacb_mat(nrows(a), ncols(b))
  ccall((:acb_mat_mul, libarb), Cint,
        (Ref{nacb_mat}, Ref{nacb_mat}, Ref{nacb_mat}, Int),
        z, a, b, p)
  return z
end

function Base.hcat(a::nacb_mat, b::nacb_mat)
  @assert nrows(a) == nrows(b)
  z = nacb_mat(nrows(a), ncols(a) + ncols(b))
  t = nacb()
  for i in 1:nrows(a)
    for j in 1:ncols(a)
      getindex!(t, a, i, j)
      setindex!(z, t, i, j)
    end
    for j in 1:ncols(b)
      getindex!(t, b, i, j)
      setindex!(z, t, i, ncols(a) + j)
    end
  end
  return z
end

function overlaps(a::narb, b::narb)
  r = ccall((:arb_overlaps, libarb), Cint,
            (Ref{narb}, Ref{narb}),
            a, b)
  return r != 0
end

function overlaps(a::nacb, b::nacb)
  r = ccall((:acb_overlaps, libarb), Cint,
            (Ref{nacb}, Ref{nacb}),
            a, b)
  return r != 0
end

function swap!(a::nacb, b::nacb)
  ccall((:acb_swap, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}),
        a, b)
end

function inv_series(a::narb_poly, ord::Int, p::Int)
  z = narb_poly()
  ccall((:arb_poly_inv_series, libarb), Nothing,
        (Ref{narb_poly}, Ref{narb_poly}, Int, Int),
        z, a, ord, p)
  return z  
end

function log_series(a::narb_poly, ord::Int, p::Int)
  z = narb_poly()
  ccall((:arb_poly_log_series, libarb), Nothing,
        (Ref{narb_poly}, Ref{narb_poly}, Int, Int),
        z, a, ord, p)
  return z  
end

function log1p_series(a::narb_poly, ord::Int, p::Int)
  z = narb_poly()
  ccall((:arb_poly_log1p_series, libarb), Nothing,
        (Ref{narb_poly}, Ref{narb_poly}, Int, Int),
        z, a, ord, p)
  return z  
end

function exp_series(a::narb_poly, ord::Int, p::Int)
  z = narb_poly()
  ccall((:arb_poly_exp_series, libarb), Nothing,
        (Ref{narb_poly}, Ref{narb_poly}, Int, Int),
        z, a, ord, p)
  return z
end

function pow_series(a::narb_poly, b::narb, ord::Int, p::Int)
  z = narb_poly()
  ccall((:arb_poly_pow_arb_series, libarb), Nothing,
        (Ref{narb_poly}, Ref{narb_poly}, Ref{narb}, Int, Int),
        z, a, b, ord, p)
  return z  
end

function mul!(z::narb_poly, a::narb_poly, b::narb, p::Int)
  ccall((:arb_poly_scalar_mul, libarb), Nothing,
        (Ref{narb_poly}, Ref{narb_poly}, Ref{narb}, Int),
        z, a, b, p)
  return z  
end

function mul(a::narb_poly, b::narb, p::Int)
  return mul!(narb_poly(), a, b, p)
end

function div!(z::narb_poly, a::narb_poly, b::narb, p::Int)
  ccall((:arb_poly_scalar_div, libarb), Nothing,
        (Ref{narb_poly}, Ref{narb_poly}, Ref{narb}, Int),
        z, a, b, p)
  return z  
end



function Base.factorial(a::Int, p::Int)
  z = narb()
  ccall((:arb_fac_ui, libarb), Nothing,
        (Ref{narb}, UInt, Int),
        z, UInt(a), p)
  return z
end

function add_error!(z::nacb, a::narb)
  GC.@preserve z begin
    t = ccall((:acb_real_ptr, libarb), Ptr{narb}, (Ref{nacb}, ), z)
    ccall((:arb_add_error, libarb), Nothing, (Ptr{narb}, Ref{narb}), t, a)
    t = ccall((:acb_imag_ptr, libarb), Ptr{narb}, (Ref{nacb}, ), z)
    ccall((:arb_add_error, libarb), Nothing, (Ptr{narb}, Ref{narb}), t, a)
  end
  return z
end

nrows(a::fmpz_mat) = a.r

ncols(a::fmpz_mat) = a.c


function zero!(z::fmpz_mat)
  ccall((:fmpz_mat_zero, libflint), Nothing,
        (Ref{fmpz_mat}, ),
        z)
end

function one!(z::fmpz_mat)
  ccall((:fmpz_mat_one, libflint), Nothing,
        (Ref{fmpz_mat}, ),
        z)
end

function one!(z::fmpz)
  ccall((:fmpz_set_ui, libflint), Nothing,
        (Ref{fmpz}, UInt),
        z, 1)
end


function magnitude(a::nacb)
  z = nmag()
  ccall((:acb_get_mag, libarb), Nothing,
        (Ref{nmag}, Ref{nacb}),
        z, a)
  return z
end

function isspecial(a::nmag)
  return 0 != ccall((:mag_is_special, libarb), Cint,
                    (Ref{nmag},),
                    a)
end

function iszero(a::nmag)
  return 0 != ccall((:mag_is_zero, libarb), Cint,
                    (Ref{nmag},),
                    a)
end

function Base.min(a::nmag, b::nmag)
  t = ccall((:mag_cmp, libarb), Cint,
            (Ref{nmag}, Ref{nmag}),
            a, b)
  return t < 0 ? a : b
end

function Base.div(a::nmag, b::nmag)
  z = nmag()
  ccall((:mag_div, libarb), Nothing,
        (Ref{nmag}, Ref{nmag}, Ref{nmag}),
        z, a, b)
  return z
end

function Base.exponent(a::nmag)
  z = fmpz()
  # dirty
  ccall((:fmpz_set, :libflint), Nothing,
        (Ref{fmpz}, Ref{nmag}),
        z, a)
  return z
end

#### special functions required for field F with elem_type T ###################

function numerator_ring(::FlintRationalField)
  return FlintIntegerRing()
end

# integer check
function isinteger_with_integer(x::fmpq)
  if isone(denominator(x))
    return true, numerator(x)
  else
    return false, ZZ()
  end
end

# ball approximation
function narb(x::fmpq, p::Int)
  z = narb()
  ccall((:arb_set_fmpq, libarb), Nothing,
        (Ref{narb}, Ref{fmpq}, Int),
        z, x, p)
  return z
end

function nacb(x::fmpq, p::Int)
  z = nacb()
  ccall((:acb_set_fmpq, libarb), Nothing,
        (Ref{nacb}, Ref{fmpq}, Int),
        z, x, p)
  return z
end

function sin_pi(a::fmpq, p::Int)
  z = narb()
  ccall((:arb_sin_pi_fmpq, libarb), Nothing,
        (Ref{narb}, Ref{fmpq}, Int),
        z, a, p)
  return z
end

function gamma(a::nacb, p::Int)
  z = nacb()
  ccall((:acb_gamma, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Int),
        z, a, p)
  return z
end

function rgamma(a::nacb, p::Int)
  z = nacb()
  ccall((:acb_rgamma, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Int),
        z, a, p)
  return z
end

function gamma(a::fmpq, p::Int)
  z = narb()
  ccall((:arb_gamma_fmpq, libarb), Nothing,
        (Ref{narb}, Ref{fmpq}, Int),
        z, a, p)
  return z
end

function rgamma(a::fmpq, p::Int)
  return inv(gamma(a,p),p)
end

function rising_factorial(a::fmpq, n::Union{fmpz, Int}, p::Int)
  z = narb()
  ccall((:arb_rising_fmpq_ui, libarb), Nothing,
        (Ref{narb}, Ref{fmpq}, UInt, Int),
        z, a, UInt(n), p)
  return z
end

function hypergeometric_2f1_regularized(a, b, c, x, Φ::Int)
  A = nacb(a, Φ)
  B = nacb(b, Φ)
  C = nacb(c, Φ)
  X = nacb(x, Φ)
  z = nacb()
  ccall((:acb_hypgeom_2f1, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Ref{nacb}, Ref{nacb}, Ref{nacb}, Cint, Int),
        z, A, B, C, X, 1, Φ)
  return z
end

function addmul!(z::fmpz, x::fmpz, y::fmpz)
   ccall((:fmpz_addmul, libflint), Nothing,
         (Ref{fmpz}, Ref{fmpz}, Ref{fmpz}),
         z, x, y)
   return z
end

function addmul!(z::fmpq, x::fmpq, y::fmpq)
   ccall((:fmpq_addmul, libflint), Nothing,
         (Ref{fmpq}, Ref{fmpq}, Ref{fmpq}),
         z, x, y)
   return z
end

function nbits(a::Int)
  return nbits(fmpz(a))
end

function nbits(a::fmpq)
  return ccall((:fmpq_height_bits, libflint), UInt,
               (Ref{fmpq},),
               a)
end

function nbits(a::fmpq_poly)
  den = fmpz_mock(a.den)
  den_bits = ccall((:fmpz_bits, libflint), UInt,
                   (Ref{fmpz_mock},),
                   den)
  num = fmpz_poly_mock(a.coeffs, a.alloc, a.length)
  num_bits = ccall((:fmpz_poly_max_bits, libflint), Int,
                   (Ref{fmpz_poly_mock},),
                   num)
  return max(den_bits, abs(num_bits)%UInt)
end

# hack
function derivative(f::FracElem)
  n = numerator(f)
  d = denominator(f)
  return (derivative(n)*d - n*derivative(d))//d^2
end


#### recurrence relations and diff equations ###################################

# equation c maintains θ on right as c0(x) + c1(x)*θ + c2(x)*θ^2 + ...
# do c = fmul*θ*fmul^-1*c + d*c where logdm = fmul'/fmul
function equ_θ_mul_add!(c::Vector, logdm, arg, d)
  FFx = parent(arg)   # F(x)
  x = FFx(gen(base_ring(FFx)))
  n = length(c)
  push!(c, zero(FFx))
  v = divexact(arg, derivative(arg))
  for i in n:-1:0
    t = (i == 0) ? FFx(0) : divexact(c[1+i-1], x)
    t += derivative(c[1+i]) - c[1 + i]*logdm
    c[1+i] = c[1+i]*d + v*t
  end
end

# equation c maintains θ on left as c0(x) + θ*c1(x) + θ^2*c2(x) + ...
# do c = c*θ + d
function equ_mul_θ_add!(c::Vector, d, Fx)
  n = length(c)
  push!(c, Fx())
  for i in n:-1:0
    c[1+i] = (i == n) ? Fx(0) : -shift_left(derivative(c[1+i]), 1)
    c[1+i] += (i == 0) ? d : c[1+i-1]
  end
end

# return either c1 - c2 or c2 - c1, doesn't matter which
function equ_sub!(c1::Vector, c2::Vector)
  if length(c1) < length(c2)
    (c1, c2) = (c2, c1)
  end
  for i in 1:length(c2)
    c1[i] -= c2[i]
  end
  return c1
end

function nacb_poly(a::fmpq_poly, Φ::Int)
  z = nacb_poly()
  ccall((:acb_poly_set_fmpq_poly, libarb), Nothing,
        (Ref{nacb_poly}, Ref{fmpq_poly}, Int),
        z, a, Φ)
  return z
end

mutable struct BiElem{T}
  isexact::Bool
  haveapprox::Bool
  exact::T
  approx::nacb
end

# TODO get rid of BiFrac
mutable struct BiFrac{T}
  isexact::Bool
  haveapprox::Bool
  num::T
  den::T
  approx::nacb
end

mutable struct BiMatrix{T}
  isexact::Bool
  haveapprox::Bool
  num::Array{T, 3}
  den::T
  approx::Array{nacb, 3}
end

function BiElem{T}(a::nacb, F) where T
  @assert T == elem_type(F)
  return BiElem{T}(false, true, zero(F), a)
end

function BiElem{T}(a::T, F) where T
  @assert T == elem_type(F)
  return BiElem{T}(true, false, a, nacb())
end

function BiFrac{T}(a::nacb, F) where T
  @assert T == elem_type(F)
  return BiFrac{T}(false, true, zero(F), zero(F), a)
end

function BiFrac{T}(n::T, d::T, F) where T
  @assert T == elem_type(F)
  return BiFrac{T}(true, false, n, d, nacb())
end

function BiMatrix{T}(m::Int, n::Int, τ::Int, F) where {T}
  @assert T == elem_type(F)
  return BiMatrix(true, false,
            [zero(F) for i in 1:m, j in 1:n, k in 1:τ], one(F),
            [nacb()  for i in 1:m, j in 1:n, k in 1:τ])
end

function get_approx(z::BiElem, Φ::Int)
  if !z.haveapprox
    z.haveapprox = true
    z.approx = nacb(z.exact, Φ)
  end
  return z.approx
end

function get_approx(z::BiFrac, Φ::Int)
  if !z.haveapprox
    z.haveapprox = true
    z.approx = div(nacb(z.num, Φ), nacb(z.den, Φ), Φ)
  end
  return z.approx
end

function force_approx(z::BiFrac, Φ::Int)
  z.isexact || return
  z.haveapprox = true
  z.isexact = false
  z.approx = div(nacb(z.num, Φ), nacb(z.den, Φ), Φ)
end

function get_approx(z::BiMatrix, Φ::Int)
  if !z.haveapprox
    z.haveapprox = true
    (m, n, τ) = size(z.num)
    den = nacb(z.den, Φ)
    for i in 1:m, j in 1:n, k in 1:τ
      z.approx[i,j,k] = div(nacb(z.num[i,j,k], Φ), den, Φ)
    end
  end
  return z.approx
end

function force_approx(z::BiMatrix, Φ::Int)
  z.isexact || return
  z.isexact = false
  get_approx(z, Φ)
end


function limit_prec(r::BiMatrix, Φ::Int)
  r.isexact || return
  if nbits(r.den) > Φ
    force_approx(r, Φ)
    return
  end
  (m, n, τ) = size(r.num)
  for i in 1:m, j in 1:n, k in 1:τ
    if nbits(r.num[i,j,k]) > Φ
      force_approx(r, Φ)
      return
    end
  end
end

function mullow(a::BiMatrix{T}, b::BiMatrix{T}, Φ::Int) where {T}
  (m, p, τ) = size(a.num)
  (p2, n, τ2) = size(b.num)
  @assert p == p2 && τ == τ2

  F = parent(a.num[1,1,1])
  znum = T[zero(F) for i in 1:m, j in 1:n, k in 1:τ]
  zapprox = nacb[zero(nacb) for i in 1:m, j in 1:n, k in 1:τ]

  if a.isexact && b.isexact
    for i in 1:m, j in 1:n, h in 1:p, k in 0:τ-1, l in 0:k
      addmul!(znum[i,j,1+k], a.num[i,h,1+l], b.num[h,j,1+k-l])
    end
    return BiMatrix{T}(true, false, znum, a.den*b.den, zapprox)
  end
  
  if !a.isexact && !b.isexact
    mullow!(zapprox, a.approx, b.approx, Φ)
  else
    t = nacb()
    if !a.isexact && b.isexact
      for i in 1:m, j in 1:n, h in 1:p, k in 0:τ-1, l in 0:k
        mul!(t, a.approx[i,h,1+l], b.num[h,j,1+k-l], Φ)
        add!(zapprox[i,j,1+k], zapprox[i,j,1+k], t, Φ)
      end
      t = nacb(b.den, Φ)
    else
      for i in 1:m, j in 1:n, h in 1:p, k in 0:τ-1, l in 0:k
        mul!(t, b.approx[h,j,1+k-l], a.num[i,h,1+l], Φ)
        add!(zapprox[i,j,1+k], zapprox[i,j,1+k], t, Φ)
      end
      t = nacb(a.den, Φ)
    end
    for i in 1:m, j in 1:n, h in 1:τ
      div!(zapprox[i,j,h], zapprox[i,j,h], t, Φ)
    end
  end
  return BiMatrix{T}(false, true, znum, zero(F), zapprox)
end

function mullow!(z::Array{nacb, 3}, a::Array{nacb, 3}, b::Array{nacb, 3}, Φ)
  (m, p, τ) = size(a)
  (p2, n, τ2) = size(b)
  @assert p == p2 && τ == τ2
  t = nacb()
  for i in 1:m, j in 1:n, h in 1:p, k in 0:τ-1, l in 0:k
    mul!(t, a[i,h,1+l], b[h,j,1+k-l], Φ)
    add!(z[i,j,1+k], z[i,j,1+k], t, Φ)
  end
  return z
end

function mullow(a::Array{nacb, 3}, b::Array{nacb, 3}, Φ)
  (m, p, τ) = size(a)
  (p2, n, τ2) = size(b)
  @assert p == p2 && τ == τ2
  return mullow!([zero(nacb) for i in 1:m, j in 1:n, k in 1:τ], a, b, Φ)
end

function mul!(z::Array{nacb, 3}, a::nacb, Φ::Int)
  for i in z
    mul!(i, i, a, Φ)
  end
end

function add!(z::Array{nacb, 3}, a::Array{nacb, 3}, b::Array{nacb, 3}, Φ::Int)
  (m, n, τ) = size(a)
  for i in 1:m, j in 1:n, k in 1:τ
    add!(z[i,j,k], a[i,j,k], b[i,j,k], Φ)
  end
  return z
end

function add(a::Array{nacb, 3}, b::Array{nacb, 3}, Φ::Int)
  (m, n, τ) = size(a)
  return nacb[add(a[i,j,k], b[i,j,k], Φ) for i in 1:m, j in 1:n, k in 1:τ]
end

function mul!(z::Array{T, 3}, a::Array{T, 3}, b::T) where T
  (m, n, τ) = size(a)
  for i in 1:m, j in 1:n, k in 1:τ
    mul!(z[i,j,k], a[i,j,k], b)
  end
end

function addmul!(z::Array{T, 3}, a::Array{T, 3}, b::T) where T
  (m, n, τ) = size(a)
  for i in 1:m, j in 1:n, k in 1:τ
    addmul!(z[i,j,k], a[i,j,k], b)
  end
end


function normalize_content(num::Matrix{PT}, den::PT) where PT
  (s, τ) = size(num)

  ZF = FlintIntegerRing()
  Zθ,_ = PolynomialRing(ZZ, "θ")


  g = content(den)
  for i in 1:s, j in 1:τ
    g = gcd(g, content(num[i,j]))
  end

  numZ = elem_type(Zθ)[map_coeffs(c-> numerator(c), divexact(num[i,j], g), parent=Zθ)
                       for i in 1:s, j in 1:τ]

  denZ::elem_type(Zθ) = map_coeffs(c-> numerator(c), divexact(den, g), parent=Zθ)

  return (numZ, denZ)
end



mutable struct SumBsCtx{ZT, PZT}
  N0::Int
  N1::Int
  eqnum::Matrix{PZT}
  eqden::PZT
  z::BiFrac{ZT}
  λnum::ZT
  λden::ZT
  Σcurr::BiMatrix{ZT}
  Σprev::BiMatrix{ZT}
  Δ::BiMatrix{ZT}
  M::BiMatrix{ZT}
end

function sum_bs_start(
  eq::Matrix{PT}, den::PT,
  z::BiElem{T},
  λ::T,
  δ::Int,
  N1::Int, N0::Int,
  Φ::Int) where {T, PT}

  F = parent(λ)
  ZF = numerator_ring(F)
  ZFθ, _ = PolynomialRing(ZF, "θ", cached=true)

  ZT = elem_type(ZF)
  PZT = elem_type(ZFθ)

  (s, τ) = size(eq)

  g = content(den)
  for i in 1:s, j in 1:τ
    g = gcd(g, content(eq[i,j]))
  end

  eqnum = PZT[map_coeffs(numerator, divexact(eq[i,j], g), parent=ZFθ) for i in 1:s, j in 1:τ]
  eqden::PZT = map_coeffs(numerator, divexact(den, g), parent=ZFθ)

  λnum::ZT = numerator(λ)
  λden::ZT = denominator(λ)

  biz = z.isexact ? BiFrac{ZT}(numerator(z.exact), denominator(z.exact), ZF) :
                    BiFrac{ZT}(z.approx, ZF)
  (Σ, M) = sum_bs_der_recursive(eqnum, eqden, biz, λnum, λden, δ, N1, N0, Φ)
  Σp = BiMatrix{ZT}(δ, s, τ, ZF)
  Δ = BiMatrix{ZT}(δ, s, τ, ZF)
  return SumBsCtx{ZT, PZT}(N0, N1, eqnum, eqden, biz, λnum, λden, Σ, Σp, Δ, M)

end


function sum_bs_continue(
  S::SumBsCtx{ZT, PZT},
  N::Int,
  Φ::Int) where {ZT, PZT}

  (δ, s, τ) = size(S.Σcurr.approx)
  (Σ, S.M, S.Δ) = sum_bs_der_continue(S.eqnum, S.eqden, S.z, S.λnum, S.λden,
                                      δ, N, S.N1, S.N0, S.Σcurr, S.M, Φ)
  (S.Σcurr, S.Σprev) = (Σ, S.Σcurr)
  S.N1 = N

  return overlap_dist(get_approx(S.Σcurr, Φ),
                      get_approx(S.Δ, Φ),
                      get_approx(S.Σprev, Φ), Φ)
end

function get_sum(S::SumBsCtx, Φ::Int)
  return get_approx(S.Σcurr, Φ)
end

# let f(z) = Σ_{N0≤i<N1, 0≤j<τ} u[i,j]*z^(λ+i)*log(z)^j/j!
# return [z^(d-N0-λ)*f^(d)(z)]_{0≤d<δ} as a δxs matrix of acb_poly in Λ
function sum_bs_der(
  eq::Matrix{PT}, den::PT,  # size sxτ
  z::BiFrac{T},
  λnum::T, λden::T,
  δ::Int,
  N1::Int, N0::Int,
  Φ::Int) where {PT, T}

  return sum_bs_der_recursive(eq, den, z, λ, δ, N1, N0, Φ)
end

function sum_bs_der_continue(
  eq::Matrix{PT}, den::PT,  # size sxτ
  z::BiFrac{T},
  λnum::T, λden::T,
  δ::Int,
  hi::Int, mid::Int, lo::Int,
  Σ0::BiMatrix{T}, M0::BiMatrix{T},
  Φ::Int) where {T, PT}

  (s, τ) = size(eq)
  F = parent(λnum)
  Σ1, M1 = sum_bs_der_recursive(eq, den, z, λnum, λden, δ, hi, mid, Φ)
  Δ = mullow(Σ1, M0, Φ)
  limit_prec(Δ, Φ)
  Σ = BiMatrix{T}(δ, s, τ, F)
  if z.isexact && Δ.isexact && Σ0.isexact
    # Σ0n/Σ0d + Δn/Δd * zn^k/zd^k = Σ0n/Σ0d + Σ1/Σ1d*M0/M0d * zn^k/zd^k
    # M0d should "divide" Σ0d
    zdenk = z.den^(mid-lo)
    znumk = z.num^(mid-lo)
    Σ.den = zdenk*Σ0.den*Σ1.den
    mul!(Σ.num, Δ.num, divexact(Σ0.den, M0.den)*znumk)
    addmul!(Σ.num, Σ0.num, Σ1.den*zdenk)
    limit_prec(Σ, Φ)
    Δ.den = Δ.den*zdenk
    mul!(Δ.num, Δ.num, znumk)
  else
    force_approx(Δ, Φ)
    mul!(Δ.approx, pow(get_approx(z, Φ), mid-lo, Φ), Φ)
    Σ.isexact = false
    Σ.haveapprox = true
    add!(Σ.approx, get_approx(Σ0, Φ), Δ.approx, Φ)
  end
  M = mullow(M1, M0, Φ)
  limit_prec(M, Φ)
  return (Σ, M, Δ)
end

function sum_bs_der_recursive(
  eq::Matrix{PT}, den::PT,  # size sxτ
  z::BiFrac{T},
  λnum::T, λden::T,
  δ::Int,
  hi::Int, lo::Int,
  Φ::Int) where {T, PT}

  @assert hi>lo

  if hi > lo + (z.isexact ? 10 : 1)
    mid = cld(hi+lo,2)
    (Σ0, M0) = sum_bs_der_recursive(eq, den, z, λnum, λden, δ, mid, lo, Φ)
    (Σ, M, Δ) = sum_bs_der_continue(eq, den, z, λnum, λden, δ, hi, mid, lo, Σ0, M0, Φ)
    return (Σ, M)
  end

  (s, τ) = size(eq)
  F = parent(λnum)

  Σ = BiMatrix{T}(δ, s, τ, F)
  M = BiMatrix{T}(s, s, τ, F)

  de = evaluate(den, lo)
  M.den = de
  for j in 1:s, k in 1:τ
    M.num[1,j,k] = evaluate(eq[j,k], lo)
  end
  for i in 2:s
    M.num[i,i-1,1] = deepcopy(de)
  end

  ff = T[λden^(δ-1)]  # λden^(δ-1)*(Λ+λ+n-0)*(Λ+λ+n-1)*...*(Λ+λ+n-(d-1))
  mul!(Σ.den, ff[1+0], de)
  for j in 1:s, k in 1:τ
    mul!(Σ.num[1,j,k], ff[1+0], M.num[1,j,k])
  end
  for d in 1:δ-1
    shift = (lo-(d-1))*λden+λnum
    # ff *= lo-(d-1)+λ + Λ
    for i in 1:length(ff)
      divexact!(ff[i], ff[i], λden)
    end
    for i in length(ff):-1:1
      mul!(ff[i], ff[i], shift)
      i == 1 || addmul!(ff[i], ff[i-1], λden)
    end
    length(ff) >= τ || push!(ff, λden^(δ-1))
    for j in 1:s, k in 0:τ-1, l in 0:k
      addmul!(Σ.num[1+d,j,1+k], ff[1+l], M.num[1,j,1+k-l])
    end
  end

  if lo+1==hi
    return (Σ, M)
  end

  @assert z.isexact

  tt = [zero(F) for k in 1:τ]
  for n in lo+1:hi-1
    # update M
    de = evaluate(den, n)
    ne = T[evaluate(eq[i,k], n) for i in 1:s, k in 1:τ]
    for j in 1:s
      # update column j
      for k in 0:τ-1
        zero!(tt[1+k])
      end
      for k in 0:τ-1, h in 0:k
        addmul!(tt[1+k], ne[1,1+h], M.num[1,j,1+k-h])
      end
      for i in s:-1:2
        for k in 0:τ-1, h in 0:k
          addmul!(tt[1+k], ne[i,1+h], M.num[i,j,1+k-h])
        end
        for k in 0:τ-1
          mul!(M.num[i,j,1+k], de, M.num[i-1,j,1+k])
        end
      end
      for k in 0:τ-1
        (M.num[1,j,1+k], tt[1+k]) = (tt[1+k], M.num[1,j,1+k])
      end
    end
    mul!(M.den, M.den, de)

    # update Σ
    mul!(Σ.den, Σ.den, z.den*de)
    mul!(Σ.num, Σ.num, z.den*de)
    # Σ.num += zn^(n-lo)*ff*M.num
    ff = T[λden^(δ-1)*z.num^(n-lo)]  # zn^(n-lo)*λd^(δ-1)*(Λ+λ+n-0)*(Λ+λ+n-1)*...*(Λ+λ+n-(d-1))
    for j in 1:s, k in 1:τ
      addmul!(Σ.num[1,j,k], ff[1+0], M.num[1,j,k])
    end
    for d in 1:δ-1
      shift = (n-(d-1))*λden + λnum
      # ff *= n-(d-1)+λ + Λ
      for i in 1:length(ff)
        divexact!(ff[i], ff[i], λden)
      end
      for i in length(ff):-1:1
        mul!(ff[i], ff[i], shift)
        i == 1 || addmul!(ff[i], ff[i-1], λden)
      end
      length(ff) >= τ || push!(ff, λden^(δ-1)*z.num^(n-lo))
      for j in 1:s, k in 0:τ-1, l in 0:k
        addmul!(Σ.num[1+d,j,1+k], ff[1+l], M.num[1,j,1+k-l])
      end
    end
  end

  return (Σ, M)
end


mutable struct bimatrix
  isexact::Bool
  p::fmpz_mat
  p_approx::narb_mat
  q::fmpz_mat
  q_approx::narb_mat
  r::fmpz
end

function bimatrix(d::Int)
  return bimatrix(true,
                  zero_matrix(ZZ, d, d),
                  narb_mat(d, d),
                  zero_matrix(ZZ, d, d),
                  narb_mat(d, d),
                  zero(ZZ))
end

function force_approx!(z::bimatrix, p::Int)
  z.isexact || return
  z.isexact = false
  d = nrows(z.q)
  for i in 1:d, j in 1:d
    z.p_approx[i,j] = narb(z.p[i,j]//z.r, p)
    z.q_approx[i,j] = narb(z.q[i,j]//z.r, p)
  end
end

function limit_prec!(z::bimatrix, p::Int)
  z.isexact || return
  p += 20
  if nbits(z.r) > p
    force_approx!(z, p)
  end
  d = nrows(z.q)
  for i in 1:d, j in 1:d
    if nbits(z.p[i,j]) > p || nbits(z.q[i,j]) > p
      force_approx!(z, p)
      return
    end
  end
end


function sum_bs_recursive!(
  z::bimatrix,
  k1::Int, k2::Int,
  num::Vector{fmpz_poly}, den::fmpz_poly,
  prec::Int)

  d = length(num)

  if k2-k1<=10
    z.isexact = true
    zero!(z.p)
    one!(z.q)
    one!(z.r)
    for k in k1:k2-1
      q2 = evaluate(den, k)
      if d>1
        e = [evaluate(num[i], k) for i in 1:d]
        for j in 1:d
          tt = q2*z.q[d,j]
          z.q[d,j] = e[1]*z.q[1,j] + e[d]*z.q[d,j]
          for i in 2:d-1
            z.q[d,j] += e[i]*z.q[i,j]
            z.q[i-1,j] = q2*z.q[i,j]
          end
          z.q[d-1,j] = tt
        end
      else
        z.q[1,1] *= evaluate(num[1], k)
      end
      mul!(z.p, q2)
      add!(z.p, z.p, z.q)
      z.r *= q2
    end
    g = z.r
    for i in 1:d, j in 1:d
      g = gcd(g, z.q[i,j])
      g = gcd(g, z.p[i,j])
    end
    z.p = divexact(z.p, g)
    z.q = divexact(z.q, g)
    z.r = divexact(z.r, g)
  else
    k12 = div(k1+k2,2)
    z1 = bimatrix(d)
    sum_bs_recursive!(z1, k1, k12, num, den, prec)
    z2 = bimatrix(d)
    sum_bs_recursive!(z2, k12, k2, num, den, prec)
    limit_prec!(z1, prec)
    limit_prec!(z2, prec)
    z.isexact = z1.isexact && z2.isexact
    if z.isexact
      mul!(z.p, z2.p, z1.q)
      addmul!(z.p, z1.p, z2.r)
      mul!(z.q, z2.q, z1.q)
      mul!(z.r, z1.r, z2.r)
    else
      force_approx!(z1, prec)
      force_approx!(z2, prec)
      mul!(z.p_approx, z2.p_approx, z1.q_approx, prec)
      add!(z.p_approx, z.p_approx, z1.p_approx, prec)
      mul!(z.q_approx, z2.q_approx, z1.q_approx, prec)
    end
  end
end


function rescale(eq::Vector{fmpq_poly})
  g = content(eq[1])
  for i in 2:length(eq)
    g = gcd(g, content(eq[i]))
  end
  Zx,x = PolynomialRing(ZZ, "x")
  return elem_type(Zx)[ map_coeffs(divexact(eq[i], g), parent=Zx) do c
                          @assert isone(denominator(c))
                          numerator(c)
                        end
                      for i in 1:length(eq)]  
end

# TODO get rid of this
function sum_bs(
  k::Int,
  eq::Vector{fmpq_poly},
  iv::Vector,
  p::Int)

  d = length(eq) - 1
  req = rescale(eq)

  s = zero(narb)
  if k > d
    z = bimatrix(d)
    sum_bs_recursive!(z, d, k, req[d+1:-1:2], -req[1], p)
    force_approx!(z, p)
    t = zero(narb)
    for i in 1:d
      add!(s, s, mul!(t, z.p_approx[d,i], iv[i], p), p)
    end
  end

  for i in min(k,d):-1:1
    add!(s, s, iv[i], p+20)
  end
  s
end

#=
 a0(n)*A[n] + a1(n)*A[n-1] + a2(n)*A[n-2] + a3(n)*A[n-3] == 0

 b0(n)*B[n] + b1(n)*B[n-1] == 0

 with C[n] = A[n]*B[n]

  a0(n)*b0(n)*b0(n-1)*b0(n-2)*C[n]
 -a1(n)*b1(n)*b0(n-1)*b0(n-2)*C[n-1]
 +a2(n)*b1(n)*b1(n-1)*b0(n-2)*C[n-2]
 -a3(n)*b1(n)*b1(n-1)*b1(n-2)*C[n-3] == 0
=#
function hadamard_product1(A::Vector{S}, B0::S, B1::S) where S
  Fθ = parent(B0)
  θ = gen(Fθ)
  Ad = length(A)-1
  C = elem_type(Fθ)[]
  t = prod(evaluate(B0, θ-i+1) for i in 1:Ad)
  push!(C, A[1+0]*t)
  g = C[end]
  for i in 1:Ad
    t = divexact(t, evaluate(B0, θ-i+1))
    t = -t*evaluate(B1, θ-i+1)
    push!(C, A[1+i]*t)
    g = gcd(g, C[end])
  end
  for i in 0:Ad
    C[1+i] = divexact(C[1+i], g)
  end
  return C
end

function hadamard_product(A::Vector{S}, B::Vector{S}) where S
  ra = length(A)-1
  rb = length(B)-1
  @assert ra > 0 && rb > 0
  if rb == 1
    return hadamard_product1(A, B[1], B[2])
  elseif ra == 1
    return hadamard_product1(B, A[1], A[2])
  end

  # TODO oh my gosh is this slow

  Fθ = parent(B[1])
  FFθ = FractionField(Fθ)
  θ = gen(Fθ)
  rc = ra*rb

  u = elem_type(FFθ)[FFθ(i==j ? 1 : 0) for i in 0:rc, j in 0:ra-1]
  for i in ra:rc
    f = [-evaluate(A[1+j], θ+i)//evaluate(A[1], θ+i) for j in 1:ra]
    for j in 1:ra, k in 1:ra
      u[1+i,j] += f[k]*u[i+1-k,j]
    end
  end

  v = elem_type(FFθ)[FFθ(i==j ? 1 : 0) for i in 0:rc, j in 0:rb-1]
  for i in rb:rc
    f = [-evaluate(B[1+j], θ+i)//evaluate(B[1], θ+i) for j in 1:rb]
    for j in 1:rb, k in 1:rb
      v[1+i,j] += f[k]*v[i+1-k,j]
    end
  end

  M = zero_matrix(FFθ, rc, rc+1)
  for i in 0:rc
    for j1 in 0:ra-1, j2 in 0:rb-1
      M[1+j1*rb+j2,1+i] = u[1+i,1+j1]*v[1+i,1+j2]
    end
  end

  (nullity, NS) = nullspace(M)
  @assert nullity == 1

  c = [(evaluate(numerator(NS[1+i,1]),θ-rc), evaluate(denominator(NS[1+i,1]),θ-rc)) for i in rc:-1:0]
  cont = zero(Fθ)
  den = one(Fθ)
  for ci in c
    cont = gcd(cont, ci[1])
    den = lcm(den, ci[2])
  end
  c = elem_type(Fθ)[divexact(ci[1], cont)*divexact(den, ci[2]) for ci in c]

  return c
end

function transpose_vars(Px::Vector, Fθ)
  θ = gen(Fθ)
  Pθ = elem_type(Fθ)[]
  for i in (length(Px)-1):-1:0
    for j in 0:degree(Px[1+i])
      while j >= length(Pθ)
        push!(Pθ, zero(Fθ))
      end
      setcoeff!(Pθ[1+j], i, coeff(Px[1+i], j))
    end
  end
  return Pθ
end

function cauchy_product(Aθ::Vector, Bθ::Vector)
  Fθ = parent(Bθ[1])
  θ = gen(Fθ)
  F = base_ring(Fθ)
  Fx, x = PolynomialRing(F, "x")
  Ax = transpose_vars(Aθ, Fx)
  Bx = transpose_vars(Bθ, Fx)

  @assert length(Bx) == 2 # only thing that's needed and only thing that works

  xlogdbx = -(Bx[1+0]+x*derivative(Bx[1+1]))//(Bx[1+1])

  ra = length(Ax) - 1
  FFx = parent(xlogdbx)
  F = elem_type(FFx)[FFx(Ax[1+ra])]
  for k in ra-1:-1:0
    F = elem_type(FFx)[
          ((i<length(F)) ? (x*derivative(F[1+i])-F[1+i]*xlogdbx) : zero(FFx)) +
          ((i > 0) ? F[1+i-1] : FFx(Ax[1+k]))
        for i in 0:length(F)]
  end

  # cancel content while θ is on the right
  c = map(p -> (numerator(p), denominator(p)), F)
  cont = zero(Fx)
  den = one(Fx)
  for ci in c
    cont = gcd(cont, ci[1])
    den = lcm(den, ci[2])
  end
  c = elem_type(Fx)[divexact(ci[1], cont)*divexact(den, ci[2]) for ci in c]

  # move θ to the left
  Cx = elem_type(Fx)[]
  for ci in reverse(c)
    equ_mul_θ_add!(Cx, ci, Fx)
  end
  return transpose_vars(Cx, Fθ)
end

function initial_values(A, v::Vector{T}, len::Int) where T
  iv = v # eh, clobber input
  while length(iv)<len
    n = length(iv)
    push!(iv, -sum(iv[n-i]*evaluate(A[2+i],n) for i in 0:min(length(A)-2,n-1))//evaluate(A[1],n))
  end
  return iv
end

function initial_values(A, v::Vector{T}) where T
  return initial_values(A, v, length(A)-1)
end


#### funny stuff for the evaluation of z^a*log(z)^b ############################

# currently two branches of log are supported. not much need for more
mutable struct zBranchInfo{T}
  λ::T
  biz::BiElem{T}
  invz_isset::Bool  # true: -π < im(log(z)) ≤ π,  false: -π ≤ im(log(z)) < π
  invz::nacb        # ignored unless invz_isset = true
  logz_isset::Bool
  logz::nacb
  logzpow::Vector{nacb}   # entry [j] is log(z)^j
  zlogzpow::Vector{nacb}  # entry [j] is z*log(z)^j
end

function zBranchInfo{T}(λ::T, z) where T
  return zBranchInfo{T}(λ, BiElem{T}(z, parent(λ)), false, zero(nacb), false, zero(nacb), nacb[], nacb[])
end

function zBranchInfo{T}(λ::T, z::nacb, invz::nacb) where T
  return zBranchInfo{T}(λ, BiElem{T}(z, parent(λ)), true, invz, false, zero(nacb), nacb[], nacb[])
end

function get_approx(Z::zBranchInfo, Φ::Int)
  return get_approx(Z.biz, Φ)
end

function set_logz(Z::zBranchInfo, Φ::Int)
  Z.logz_isset && return
  if Z.invz_isset
    Z.logz = neg!(log(Z.invz, Φ))
  else
    Z.logz = log(get_approx(Z, Φ), Φ)
  end
end

# z^(λ+i)
function zpow(Z::zBranchInfo, i::Int, Φ::Int)
  if Z.invz_isset
    pow(Z.invz, -Z.λ-i, Φ)
  else
    pow(get_approx(Z, Φ), Z.λ+i, Φ)
  end
end

# given a bound on |z|, compute z^a*log(z)^j/j!
# at z = 0: z^a*log^j(z) = nan if re(a) < 0
#                          0   elseif re(a) > 0
#                          1   elseif a = 0 && j == 0
#                          nan else
# TODO: get this working for T != fmpq
function fancypow_at_0(zmag::narb, a::T, j::Int, p::Int) where T <: fmpq
  r = zero(nacb)
  if a < 0
    indeterminate!(r)
  elseif a > 0
    if iszero(zmag)
      zero!(r)
    else
      s = narb()
      t = narb()
      # |z|^a * |log(|z|) +- Pi|^j/j!
      pow!(s, zmag, a, p)
      log!(t, zmag, p)
      abs!(t, t)
      add!(t, t, 4, p)
      pow!(t, t, UInt(j), p)
      div!(t, t, factorial(j, p), p)
      mul!(s, s, t, p)
      add_error!(r, s)
    end
  elseif a == 0 && j == 0
    one!(r)
  else
    indeterminate!(r)
  end
  return r
end

# z^(λ+i)*log(z)^j/j!
function zpowlogpow(Z::zBranchInfo, i::Int, j::Int, Φ::Int)
  @assert j≥0
  if j<1
    return zpow(Z, i, Φ)
  end
  set_logz(Z, Φ)
  if isfinite(Z.logz)
    return div(mul(zpow(Z, i, Φ), pow(Z.logz, j, Φ), Φ), factorial(j, Φ), Φ)
  else
    return fancypow_at_0(abs_upper(get_approx(Z, Φ), Φ), Z.λ+i, j, Φ)
  end
end

# log(z)^j
function logzpow(Z::zBranchInfo, j::Int, Φ::Int)
  @assert j≥0
  j<1 && return one(nacb)
  set_logz(Z, Φ)
  while j>length(Z.logzpow)
    i = 1+length(Z.logzpow)
    push!(Z.logzpow, pow(Z.logz, i, Φ))
  end
  return Z.logzpow[j]
end

# z*log(z)^j
function zlogzpow(Z::zBranchInfo, j::Int, Φ::Int)
  @assert j≥0
  z = get_approx(Z, Φ)
  j<1 && return z
  set_logz(Z, Φ)
  while j>length(Z.zlogzpow)
    i = 1+length(Z.zlogzpow)
    if isfinite(Z.logz)
      res = mul(z, pow(Z.logz, i, Φ), Φ)
    else
      res = fancypow_at_0(abs_upper(z, Φ), fmpq(1), j, Φ)
    end
    push!(Z.zlogzpow, res)
  end
  return Z.zlogzpow[j]
end

# sum_{j<τ} |(z+ε)^δ*log(z+ε)|/j!  mod ε^δ
function logzpoly_series(Z::zBranchInfo, δ::Int, τ::Int, Φ::Int)
  @assert δ>0
  c = Vector{fmpq}[fmpq[QQ(binomial(ZZ(δ), ZZ(k)))] for k in 0:δ-1]
  ff = [abs(mul(get_approx(Z, Φ), c[1+k][1], Φ), Φ) for k in 0:δ-1]
  jfac = fmpz(1)
  for j in 1:τ-1
    mul!(jfac, jfac, j)
    for k in δ-1:-1:0
      pushfirst!(c[1+k], zero(QQ))
      for l in 1:k, m in 0:j-1
          c[1+k][1+m] += (-1)^(l-1)*c[1+k-l][1+m]//l
      end
    end
    for k in 0:δ-1
      tt = mul(get_approx(Z, Φ), c[1+k][1+0], Φ)
      for m in 1:j
        add!(tt, tt, mul(zlogzpow(Z, m, Φ), c[1+k][1+m], Φ), Φ)
      end
      add!(ff[1+k], ff[1+k], div(abs(tt, Φ), jfac, Φ), Φ)
    end
  end
  for k in 0:δ-2
    mul!(ff[1+k], ff[1+k], pow(abs(get_approx(Z, Φ), Φ), UInt(δ-1-k), Φ), Φ)
  end
  return narb_poly(ff)
end

# z^(λ+n)*Σ_{j<τ}e_j*Λ^(τ-1-j) with Λ^(τ-1-j) replaced by log(z)^j/j!
function zpowlogpoly_eval(Z::zBranchInfo, n::Int, e::Vector{nacb}, Φ::Int)
  τ = length(e)
  if τ<2
    return mul(zpow(Z, n, Φ), e[1+0], Φ)
  end
  set_logz(Z, Φ)
  j = τ-1
  ok = isfinite(Z.logz)
  s = mul(e[1+τ-1-j], ok ? logzpow(Z, j, Φ) : zlogzpow(Z, j, Φ), Φ)
  t = nacb()
  while (j-=1)>0
    div!(t, s, j+1, Φ)
    mul!(s, e[1+τ-1-j], ok ? logzpow(Z, j, Φ) : zlogzpow(Z, j, Φ), Φ)
    add!(s, s, t, Φ)
  end
  add!(s, s, ok ? e[1+τ-1] : mul!(t, e[1+τ-1], get_approx(Z, Φ), Φ),  Φ)
  return mul(zpow(Z, n-!ok, Φ), s, Φ)
end




#=
with F[m,k] = 2f1r(1, a+m, 1+σ+a+m+k, z) we already have to use Feq for the
recursion on k. Could try using a recursion on m too:
  {F[m,0], F[m,1]} = 1/(m+a-1)*M.{F[m-1,0], F[m-1,1]}
where M = {{1, -1-σ}, {(z-1)/(z*(m+σ+a)), 1/z + (-1-σ)/(m+σ+a)}}
Since arb's 2f1r is pretty good in this case, the outlook is not so good.
=#
function heuristic_hyp(
  a::Vector{T}, b::Vector{T},
  z::T,
  m::Int, # number of terms for plain sum
  N::Int, # number of terms for tail sum
  Φ::Int) where T

  σ::T = sum(b) - sum(a)
  q = length(b)
  @assert q>1 && q+1==length(a)
  F = parent(a[end])
  Fθ, θ = PolynomialRing(F, "θ")
  Φ += 20

  Φp = Φ+10+2*nbits(m)
  S = sum_bs_start(fill(prod(θ-1+i for i in a),1,1), θ*prod(θ-1+i for i in b),
                   BiElem{T}(z, F), zero(F), 1, m, 1, Φp)
  plain = add(mul(get_sum(S, Φp)[1,1,1], z, Φp), 1, Φp+60)

  Aeq = buehring_equ(a, b, 0, θ)
  if isone(z)
    AFeq = hadamard_product(Aeq, [(θ+σ)*(θ-1+σ+a[end]+m), -(θ-1+σ)])
    s = length(AFeq)-1
    iv = initial_values(AFeq, [F(1)], s)
    Φt = Φ+4*length(AFeq)+4*nbits(N)
    S = sum_bs_start(reshape([-AFeq[1+i] for i in 1:s],s,1), AFeq[1],
                     BiElem{T}(one(F), F), zero(F), 1, N, s, Φt)
    Σa = get_sum(S, Φt)
    tail = zero(nacb)
    ivsum = zero(F)
    for i in 0:s-1
      add!(ivsum, ivsum, iv[1+s-i-1])
      add!(tail, tail, mul(Σa[1,1+i,1], iv[1+s-i-1], Φt), Φt)
    end
    add!(tail, tail, nacb(ivsum, Φt), Φt)
    mul!(tail, tail, mul(rgamma(σ+a[end]+m, Φ), inv(σ), Φ), Φ)
  else
    Feq = [z*(θ+σ)*(θ-1+σ+a[end]+m), -(θ-1+σ-(1-z)*(2*θ-2+2*σ+a[end]+m)), -Fθ(1-z)]
    AFeq = hadamard_product(Aeq, Feq)
    s = length(AFeq)-1
    Φt = Φ+4*length(AFeq)+4*nbits(N)
    iv = initial_values(Aeq, [F(1)], s)
    zz = convert(Complex{Float64}, nacb(z, 60))
    Φ2f1 = Φ + min(Φ, trunc(Int, 0.25*m*nbits(m)*abs(1-zz)))
    Fk0 = hypergeometric_2f1_regularized(F(1), a[end]+m, 1+σ+a[end]+m+0, z, Φ2f1)
    Fk1 = hypergeometric_2f1_regularized(F(1), a[end]+m, 1+σ+a[end]+m+1, z, Φ2f1)
    iv2 = nacb[mul(Fk0, iv[1+0], Φ), mul(Fk1, iv[1+1], Φ)]
    M = (F(1), F(0), F(0), F(1))
    for k in 2:length(iv)-1
      f0 = evaluate(Feq[1+0], k)
      f1 = -evaluate(Feq[1+1], k)//f0
      f2 = (1-z)//f0
      M = (f1*M[1]+f2*M[3], f1*M[2]+f2*M[4], M[1], M[2])
      Fk = add(mul(Fk1, M[1], Φ), mul(Fk0, M[2], Φ), Φ)
      push!(iv2, mul(Fk, iv[1+k], Φ))
    end
    S = sum_bs_start(reshape([-AFeq[1+i] for i in 1:s],s,1), AFeq[1],
                     BiElem{T}(one(F), F), zero(F), 1, N, s, Φt)
    Σa = get_sum(S, Φt)
    tail = zero(nacb)
    ivsum = zero(nacb)
    for i in 0:s-1
      add!(ivsum, ivsum, iv2[1+i], Φt)
      add!(tail, tail, mul(Σa[1,1+s-1-i,1], iv2[1+i], Φt), Φt)
    end
    add!(tail, tail, ivsum, Φt)
  end
  f = mul(pow(nacb(z, Φ), UInt(m), Φ), rising_factorial(a[end], m, Φ), Φ)
  for i in 1:q
    mul!(f, f, div(gamma(b[i], Φ), gamma(a[i], Φ), Φ), Φ)
  end
  return add(plain, mul(f, tail, Φ), Φ)
end


# Return equation satisfied by m(x)*pFq(arg(x)) with θ on the left.
# Allowed m's are those for which logdm(x) = m'(x)/m(x) is rational.
# Equation is returned as the dense coefficient list of powers of x.
function hyp_equ(a::Vector{T}, b::Vector{T}, logdm, arg) where T
  FFx = parent(arg)   # F(x)
  Fx = base_ring(FFx) # F[x]
  F = base_ring(Fx)   # F
  @assert FFx == parent(logdm)

  # c = equation with θ on the right
  c = elem_type(FFx)[arg]
  for ai in a
    equ_θ_mul_add!(c, logdm, arg, ai-1)
  end
  c2 = elem_type(FFx)[one(FFx)]  
  equ_θ_mul_add!(c2, logdm, arg, F(0)) # "extra" denominator param 1
  for bi in b
    equ_θ_mul_add!(c2, logdm, arg, bi-1)
  end
  c = equ_sub!(c, c2)

  # cancel content while θ is on the right
  c = map(p -> (numerator(p), denominator(p)), c)
  cont = zero(Fx)
  den = one(Fx)
  for ci in c
    cont = gcd(cont, ci[1])
    den = lcm(den, ci[2])
  end
  c = elem_type(Fx)[divexact(ci[1], cont)*divexact(den, ci[2]) for ci in c]

  # move θ to the left
  Px = elem_type(Fx)[]
  for ci in reverse(c)
    equ_mul_θ_add!(Px, ci, Fx)
  end

  # transpose so θ is inside P0(θ) + P1(θ)*x + ... + Pn(θ)*x^n
  Fθ, θ = PolynomialRing(F, "θ")
  Pθ = elem_type(Fθ)[]
  for i in (length(Px)-1):-1:0
    for j in 0:degree(Px[1+i])
      while j >= length(Pθ)
        push!(Pθ, zero(Fθ))
      end
      setcoeff!(Pθ[1+j], i, coeff(Px[1+i], j))
    end
  end

  return (Pθ, Px)
end


#### evaluation of series solutions to diff equations ##########################
#=
The setup follows Mezzarobba 2019 Truncation Bounds

With θ on the left we have a differential operator

  P(x,θ) = θ^r*Px_r + ... + θ*Px_1 + Px_0,  Px_i in F[x]
         = Pθ_s*x^s + ... + Pθ_1*x + Pθ_0,  Pθ_i in F[θ]

The radius of convergence of the solutions about 0 are related to the zeros of
Px[r]: we assume that Px_r does not have a zero at x = 0. In this case series
solution to P(z,θ)u(z) = 0 can be written down as

  u(z) = Σ_{0≤i<∞, 0≤j≤K} u_{i,j} z^(λ+i)*log^j(z)/j!

for some K, and which u_{i,j} are free and which u_{i,j} determined by pervious
terms in the solution is explained in Corollary 5.4: Only the u_{i,j} for
which j is strictly less than the multiplicy of λ+i as a root of Pθ_0 are free.
Considering all λ such that λ is a root of Pθ_0 and none of λ-1, λ-2, ... are
roots of Pθ[0] gives a total of deg(Pθ_0) linearly independent solutions.

Now consider a truncation

  u_N(z) = Σ_{0≤i<N, 0≤j≤K} u_{i,j} z^(λ+i)*log^j(z)/j!

and the normalized difference y(z) = Px_r(z)*(u_N(z) - u(z)).
This y(z) satisfies L(z,θ)y(z) = Q_0(θ)q(z) where

  1. L(x,θ) := P(z,θ)/Px_r(z) = Q_0(θ) + Q_1(θ)*z + Q_2(θ)*z^2 + ...
  2. Q_0(θ) is monic in θ with degree r (and is proportional to Pθ_0)
  3. The Q_1, Q_2, ... all have degree < r
  4. Assuming none of λ+N, ..., λ+N+s-1 are roots of Q_0(θ), the normalized
     residue q(z) is of the form

        q(z) = Σ_{0≤i<s, 0≤j≤K} q_{i,j} z^(λ+N+i)*log^j(z)/j!

Let ahat(z) be a series satisfying Proposition 5.5. Compute this by noting

    Σ_{1≤j} Q_j(θ)*x^j = P(x,θ)/Px_r(x) - Q_0(θ)
                       = P(x,θ)/Px_r(x) - P(0,θ)/Px_r(0)

If we can expand the rhs as a finite linear combination

    Σ_i f_i(θ)*(some power series in RR_+[[x]])

then it suffices to bound each of the f_j(θ)/Q_0(θ) in accordance with the lhs
of 5.9. Since we are ultimately interested in hhat(z) = exp(int_0^z ahat(z)/z dz),
it makes since to expand

    Σ_{1≤j} Q_j(θ)*x^j = Σ_i f_i(θ) x*d/dx(some power series in x*RR_+[[x]])

For all the differential equations considered here, since we have few
singularities, the "some power series in x*RR_+[x]" can be taken to be

  log(1/(1-x)),  log(1/(1-x^2)),  log((1+x)/(1-x)),
  x/(1-x)^i,     x^1/(1-x^2)^i,   x^2/(1-x^2)^i        1≤i≤~s

and our hhat(z) will take the shape (1-z)^?*(1+z)^?*exp(rational function of z).

At this point, with some minor technicalities, Algorithm 6.11 computes a ghat
such that

(*)   z^-λ*(u_N(z) - u(z)) << phat(z)*z^N*ghat(z)*hhat(z)

where phat(z) is a majorant of 1/Px_r(z). Note that we have taken out the z^N
from the ghat(z). Since the lhs is a polynomial in log(z), this must be
interpreted as saying that the rhs majorizes each of the coefficients of log(z).

Given (*) how to compute bounds on |u_N(z) - u(z)| and the derivatives
|u_N'(z) - u'(z)|, |u_N''(z) - u''(z)|, ... |u_N^(δ-1)(z) - u^(δ-1)(z)| ???

The code currently computes the power series to order O(ε^δ) of

  phat(|z|+ε)*(|z|+ε)^(λ+N)*ghat(|z|+ε)*hhat(|z|+ε)*(Σ_{j<?} log^j(|z|+ε)/j!

and does an add_error! on the coefficients of ε^d/d!, but this looks suspicious.
=#

#=
Example.

Consider computing 2F1(a,b,c,z) following the "compute_f_anywhere" path.

The equation satisfied by u(x) = (1+x)^(-2*a)*2F1(a,b,c,4*x/(1+x)^2) is

  julia> R,(a,b,c)=PolynomialRing(ZZ,["a","b","c"])

  julia> F=FractionField(R)

  julia> (a,b,c)=map(F,(a,b,c))

  julia> Fx, x = PolynomialRing(F, "x")

  julia> hyp_equ([a,b], [c], -2*Fx(a)//(1+x), 4*x//(1+x)^2)

  (AbstractAlgebra.Generic.Poly{AbstractAlgebra.Generic.Frac{fmpz_mpoly}}[
    θ^2 + (c - 1)*θ,
    (-4*b + 2*c)*θ - 4*a*b + 2*a*c + 4*b - 2*c,
    -θ^2 + (-4*a + c + 3)*θ - 4*a^2 + 2*a*c + 6*a - 2*c - 2],

  AbstractAlgebra.Generic.Poly{AbstractAlgebra.Generic.Frac{fmpz_mpoly}}[
    (-4*a^2 + 2*a*c + 6*a - 2*c - 2)*x^2 + (-4*a*b + 2*a*c + 4*b - 2*c)*x,
    (-4*a + c + 3)*x^2 + (-4*b + 2*c)*x + (c - 1),
    -x^2 + 1])

which indicates

P(x, θ) = -(θ+2*a-2)*(θ+2*a-c-1)*x^2 - 2*(2*b-c)*(θ+a-1)*x + θ*(θ+c-1)
        = θ^2*(1-x)*(1+x) + ...


Now the partial fraction decomposition after integration is

  julia> Fθ,θ=PolynomialRing(F,"θ")

  julia> Fθz,z=PolynomialRing(Fθ,"z")

  julia> partial_fractions(divexact(-(θ+2*a-2)*(θ+2*a-c-1)*z^2 -
             2*(2*b-c)*(θ+a-1)*z + θ*(θ+c-1) - θ*(θ+c-1)*(1-z)*(1+z), z), 1, 1)

  ((-2*a + c + 1)*θ - 2*a^2 + a*c + 3*a - c - 1)*log(1/(1 - z^2)) +
     ((-2*b + c)*θ - 2*a*b + a*c + 2*b - c)*log((1 + z)/(1 - z))


Since Q_0(θ) = θ*(θ+c-1), For large N, our hhat(z) is going to be about

 hhat(z) = (1/(1-z^2))^|-2*a+c+1| * ((1+z)/(1-z))^|-2*b+c|,

and phat(z) is 1/(1-z^2).


From P(x, θ), the coefficients of u(z) = Σ_n u_n*z^n satisfy

  (n+1)*(n+c)*u_{n+1} = 2*(2*b-c)*(n+a)*u_n + (n+2*a-1)*(n+2*a-c)*u_{n-1}.

For real a,b,c, a tight asymptotic bound on any solution to this is

    u_n = O(n^(|2*b-c|+2*a-c-1)),

and the majorant method via phat(z)*hhat(z) gives

    u_n = O(n^(|2*b-c|+|2*a-c-1|)).

Thus the majorant method has overestimated n^(2*a-c-1) by n^|2*a-c-1|, which
can be a serious overestimation.

=#

# for all differential equations here the int_0^z ahat(z)/z dz series can
# expressed as a finite linear combination of the following functions
# for vector entries, the part [i] gives the coefficient
mutable struct hyp_majorant{S}
  log_coeff::S              # log(1/(1-z))
  log2_coeff0::S            # log(1/(1-z^2))
  log2_coeff1::S            # log((1+z)/(1-z))
  poly_coeffs::Vector{S}    # [z^i]
  recp_coeffs::Vector{S}    # [z/(1-z)^i]
  recp2_coeffs1::Vector{S}  # [z^1/(1-z^2)^i]
  recp2_coeffs2::Vector{S}  # [z^2/(1-z^2)^i]
end

function expressify(a::hyp_majorant; context = context)
  s = Expr(:call, :+)
  push!(s.args, Expr(:call, :*, expressify(a.log_coeff), :(log(1/(1-z)))))
  push!(s.args, Expr(:call, :*, expressify(a.log2_coeff0), :(log(1/(1-z^2)))))
  push!(s.args, Expr(:call, :*, expressify(a.log2_coeff1), :(log((1+z)/(1-z)))))
  for i in 1:length(a.recp_coeffs)
    push!(s.args, Expr(:call, :*, expressify(a.recp_coeffs[i]), :(z/(1-z)^$i)))    
  end
  for i in 1:length(a.recp2_coeffs1)
    push!(s.args, Expr(:call, :*, expressify(a.recp2_coeffs1[i]), :(z/(1-z^2)^$i)))    
  end
  for i in 1:length(a.recp2_coeffs2)
    push!(s.args, Expr(:call, :*, expressify(a.recp2_coeffs1[i]), :(z^2/(1-z^2)^$i)))    
  end
  for i in 1:length(a.poly_coeffs)
    push!(s.args, Expr(:call, :*, expressify(a.poly_coeffs[i]), :(z^$i)))    
  end
  return s
end

function coeffs(a::hyp_majorant)
  return Iterators.flatten(((a.log_coeff,),
                            (a.log2_coeff0,),
                            (a.log2_coeff1,),
                            a.poly_coeffs,
                            a.recp_coeffs,
                            a.recp2_coeffs1,
                            a.recp2_coeffs2))
end

function hyp_majorant{narb}(a::hyp_majorant{S}) where S
  return hyp_majorant{narb}(
          zero(narb),
          zero(narb),
          zero(narb),
          [zero(narb) for i in 1:length(a.poly_coeffs)],
          [zero(narb) for i in 1:length(a.recp_coeffs)],
          [zero(narb) for i in 1:length(a.recp2_coeffs1)],
          [zero(narb) for i in 1:length(a.recp2_coeffs2)])
end

function Base.show(io::IO, a::hyp_majorant)
  print(io, Nemo.AbstractAlgebra.obj_to_string(a))
end

# series expansion to O(z^ord)
function series(a::hyp_majorant{narb}, ord::Int, p::Int)
  s = zero(narb_poly)
  for j in 1:ord-1
    c = zero(narb)

#  log_coeff::S              # log(1/(1-z))
    add!(c, c, mul(a.log_coeff, fmpq(1//j), p), p)

#  log2_coeff0::S            # log(1/(1-z^2))
#  log2_coeff1::S            # log((1+z)/(1-z))
    add!(c, c, mul(iseven(j) ? a.log2_coeff0 : a.log2_coeff1, fmpq(2//j), p), p)

#  poly_coeffs::Vector{S}    # [z^i]
    for i in 1:length(a.poly_coeffs)
      add!(c, c, a.poly_coeffs[i], p)
    end

#  recp_coeffs::Vector{S}    # [z/(1-z)^i]
    for i in 1:length(a.recp_coeffs)
      d = binomial(ZZ(i-1+j-1), ZZ(j-1))
      add!(c, c, mul(a.recp_coeffs[i], d, p), p)      
    end

    if isodd(j)
#  recp2_coeffs1::Vector{S}  # [z^1/(1-z^2)^i]
      for i in 1:length(a.recp2_coeffs1)
        d = binomial(ZZ(i-1+div(j-1,2)), ZZ(div(j-1,2)))
        add!(c, c, mul(a.recp2_coeffs1[i], d, p), p)
      end
    else
#  recp2_coeffs2::Vector{S}  # [z^2/(1-z^2)^i]
      for i in 1:length(a.recp2_coeffs1)
        d = binomial(ZZ(i-1+div(j-2,2)), ZZ(div(j-2,2)))
        add!(c, c, mul(a.recp2_coeffs2[i], d, p), p)
      end
    end

    setcoeff!(s, j, c)
  end
  return s
end

# evaluate it at z = z
function eval_series(a::hyp_majorant{narb}, z::narb_poly, ord::Int, p::Int)
  l1 = neg!(log1p_series(neg(z), ord, p))
  l2 = log1p_series(z, ord, p)

#  log_coeff::S              # log(1/(1-z))
  s = iszero(a.log_coeff) ? zero(narb_poly) : mul(l1, a.log_coeff, p)

#  log2_coeff0::S            # log(1/(1-z^2))   = l1 - l2
#  log2_coeff1::S            # log((1+z)/(1-z)) = l1 + l2
  iszero(a.log2_coeff1) || add!(s, s, mul(l1, add(a.log2_coeff1, a.log2_coeff0, p), p), p)
  iszero(a.log2_coeff0) || add!(s, s, mul(l2, sub(a.log2_coeff1, a.log2_coeff0, p), p), p)

#  poly_coeffs::Vector{S}    # [z^i]
  v = z
  t = one(narb_poly)
  for i in 1:length(a.poly_coeffs)
    t = mullow(t, v, ord, p)
    iszero(a.poly_coeffs[i]) || add!(s, s, mul(t, a.poly_coeffs[i], p), p)
  end

#  recp_coeffs::Vector{S}    # [z/(1-z)^i]
  v = inv_series(sub(1, z, p), ord, p)
  t = z
  for i in 1:length(a.recp_coeffs)
    t = mullow(t, v, ord, p)
    iszero(a.recp_coeffs[i]) || add!(s, s, mul(t, a.recp_coeffs[i], p), p)
  end

#  recp2_coeffs1::Vector{S}  # [z^1/(1-z^2)^i]
  v = inv_series(sub(1, mullow(z, z, ord, p), p), ord, p)
  t = z
  for i in 1:length(a.recp2_coeffs1)
    t = mullow(t, v, ord, p)
    iszero(a.recp2_coeffs1[i]) || add!(s, s, mul(t, a.recp2_coeffs1[i], p), p)
  end

#  recp2_coeffs2::Vector{S}  # [z^2/(1-z^2)^i]
  t = mullow(z, z, ord, p)
  for i in 1:length(a.recp2_coeffs1)
    t = mullow(t, v, ord, p)
    iszero(a.recp2_coeffs1[i]) || add!(s, s, mul(t, a.recp2_coeffs1[i], p), p)
  end

  return s
end

# p += a*x^i  in F[x]
function add_coeff!(p, i, a, F)
  if iszero(a)
    return
  end
  while length(p) <= i
    push!(p, zero(F))
  end
  p[1+i] += a
end

# expansion of f/((1-z)^k1*(1+z)^k2)
function partial_fractions(f, k1::Int, k2::Int)
  Fθz = parent(f)
  Fθ = base_ring(Fθz)
  F = base_ring(Fθ)
  S = elem_type(Fθ)
  z = gen(Fθz)

  r = hyp_majorant{S}(zero(Fθ), zero(Fθ), zero(Fθ), S[], S[], S[], S[])

  if degree(f) >= k1+k2
    # divrem does not work in F[θ][z] :(
    den = (1-z)^k1*(1+z)^k2
    while degree(f) >= k1+k2
      a = divexact(coeff(f, degree(f)), coeff(den, k1+k2))
      add_coeff!(r.poly_coeffs, degree(f)-(k1+k2), a, Fθ)
      f -= a*shift_left(den, degree(f)-(k1+k2))
    end
  end

  while k1>0 || k2>0
    if k1>k2
      a = F(2)^-k2*evaluate(f, Fθ(1))
      if k1 == 1
        f = divexact(-a*(1+z)^k2+f, 1-z)
        r.log_coeff += a
      else
        a = divexact(a, k1-1)
        f = divexact(-a*(1+z)^k2*(1+(k1-2)*z)+f, 1-z)
        add_coeff!(r.recp_coeffs, k1-1, a, Fθ)
      end
      k1 -= 1
	  else
      e1 = F(2)^(k2-k1)*evaluate(f, Fθ(-1))
      f = (1-z)^(k2-k1)*f
      e2 = evaluate(f, Fθ(1))
      a = divexact(e2+e1, k2 == 1 ? 4 : 4*k2-4);
      b = divexact(e2-e1, k2 == 1 ? 4 : 4*k2-4);
      if k2 == 1
        f = divexact(-2*(a+b*z)+f, 1-z^2)
        r.log2_coeff1 += a
        r.log2_coeff0 += b
      else
        f = divexact(a*(-1+(3-2*k2)*z^2)-2*b*(z+(-2+k2)*z^3)+f,1-z^2)
        add_coeff!(r.recp2_coeffs1, k2-1, a, Fθ)
        add_coeff!(r.recp2_coeffs2, k2-1, b, Fθ)
      end
      k1 = k2 = k2 - 1
    end
  end

  return r
end

# Algorithm 7.1
# bound Σ_{j<τ} |n*[θ^j]f(n+θ)/g(n+θ)| for n0≤n≤n1
# it is good if g does not vanish at any integer in [n0, n1]
# pass n1 < 0 for n1 = ∞ :(
# should have g(θ) = prod_α (θ-α)
function fraction_bound_normal(
  f, g,
  αs::Vector{T},
  τ::Int,
  n0::fmpz, n1::fmpz,
  p::Int) where T

  d = length(αs)
  @assert d == degree(g)

  if iszero(f)
    return zero(narb)
  end

  pm_one = narb()
  ccall((:arb_zero_pm_one, libarb), Nothing,
        (Ref{narb},),
        pm_one)

  # x = interval containing [1/n1, 1/n0]
  x = nacb(mul(pm_one, 1//n0, p))
#=
  let a[b] denote coeff(a, b). Since the case of f[d-1] ≠ 0 is important,

  n*f(n+ε)/g(n+ε)

  = (n+ε)*f(n+ε)/g(n+ε) - ε*f(n+ε)/g(n+ε)

             (n+ε)*f(n+ε) - f[d-1]*g(n+ε) - ε*f(n+ε)
  = f[d-1] + ---------------------------------------
                              g(n+ε)

             Σ_{0≤i<d} (f[i-1] - f[d-i]*g[i] - ε*f[i])*(n+ε)^i
  = f[d-1] + ---------------------------------------------------
                             prod_α (n+ε-α)
  with x = 1/n
               Σ_{0≤i<d} (f[i-1] - f[d-i]*g[i] - ε*f[i])*x^(d-1-i)*(1+x*ε)^i
  = f[d-1] + x*---------------------------------------------------------------
                                      prod_α (1-x*α+x*ε)
=#
  num = zero(nacb_poly)
  t = one(nacb)
  t2 = nacb_poly(one(nacb), x)
  c0 = nacb()
  c1 = nacb()
  for i in d-1:-1:0
    mul!(c0, t, (i>0 ? coeff(f,i-1) : 0) - coeff(f,d-1)*coeff(g,i), p)
    mul!(c1, t, -coeff(f, i), p)
    add!(num, nacb_poly(c0, c1), mullow(t2, num, τ, p), p)
    mul!(t, t, x, p)
  end

  den = one(nacb_poly)
  c = one(narb)
  for α in αs
    mul!(c0, x, α, p)
    sub!(c0, 1, c0, p)
    den = mullow(nacb_poly(c0, x), den, τ, p)
    if α > 0 # TODO
      loc = numerator(ceil(α))  # TODO
      if loc <= n0
        m = abs(1-divexact(α,n0))
      elseif n1 >= 0 && n1 < loc
        m = abs(1-divexact(α,n1))
      else
        m = min(abs(1-divexact(α,loc-1)), abs(1-divexact(α,loc)))
      end
      mul!(c, c, m, p)
    end
  end
  inv!(c, c, p)

  # make c centered around 0
  mul!(c, c, pm_one, p)

  mul!(den, den, nacb(c), p)
  setcoeff!(den, 0, one(nacb))
  rat = div_series(num, den, τ, p)
  res = zero(narb)
  for i in 0:τ-1
    add!(res, res, abs(coeff(rat, i), p), p)
  end
  mul!(c, c, real(x), p)
  mul!(res, res, c, p)
  add!(res, res, abs(nacb(coeff(f,d-1), p), p), p)

  return res
end

# bound Σ_{j<τ} |n*[θ^j]f(n+θ)/(θ^-μ*g(n+θ))|
function fraction_bound_special(f, g, μ::Int, τ::Int, n::fmpz, p::Int)
  Fθ = parent(f)
  @assert Fθ == parent(g)
  θ = gen(Fθ)
  f = evaluate(f, θ+n)
  g = evaluate(g, θ+n)
  F = nacb_poly()
  G = nacb_poly()
  for j in 0:τ-1
    setcoeff!(F, j, nacb(coeff(f, j), p))
    setcoeff!(G, j, nacb(coeff(g, μ+j), p))
  end
  rat = div_series(F, G, τ, p)
  res = zero(narb)
  for j in 0:τ-1
    add!(res, res, abs(mul(coeff(rat, j), n, p), p), p)
  end
  return res  
end


function equ_bound(
  Pθ, Px,
  λ::T,
  ns::Vector{fmpz},
  αs::Vector{T},
  τ::Int,
  N::Int,
  p::Int) where T

  Fθ = parent(Pθ[1])
  θ = gen(Fθ)
  s = length(Pθ)-1
  r = length(Px)-1
  c = coeff(Px[1+r], 0)

  # Q0 = Q_0(θ) shifted by lambda
  Q0 = evaluate(divexact(Pθ[1+0],c), θ + λ)

  # f = (Pθ - Px[1+r]*Q0)/(c*x) convert to bivariate with z a.k.a x on the outside
  # f is also shifted by λ
  Fθz,z = PolynomialRing(Fθ, "z")
  f = sum(evaluate(Pθ[1+i], θ + λ)*z^i for i in 0:s)
  f = divexact(f - Px[1+r](z)*Q0, c*z)

  # also shift the roots αs of Q0
  αs = αs .- λ
  @assert Q0 == prod(θ - α for α in αs)

  # denominator should be of the form (1-x)^k1*(1+x)^k2
  x = gen(parent(Px[1+r]))
  den = divexact(Px[1+r], c)
  k1 = k2 = 0
  while degree(den) > 0
    Q, R = divrem(den, 1-x)
    if iszero(R)
      k1 += 1
      den = Q
    else
      Q, R = divrem(den, 1+x)
      if iszero(R)
        k2 += 1
        den = Q
      else
        @assert false
      end
    end
  end
  @assert isone(den)

  pf = partial_fractions(f, k1, k2)

  # start with 0's and take max's along intervals all the way from N to ∞
  maj = hyp_majorant{narb}(pf)
  curN = fmpz(N)
  curτ = τ

  if !isempty(ns)
    # collect the multiplicities of the ns
    nμ = Tuple{fmpz, Int}[]
    for n in ns
      if !isempty(nμ) && nμ[end][1] == n
        nμ[end] = (n, nμ[end][2]+1)
      else
        push!(nμ, (n, 1))
      end
    end

    for (n, μ) in nμ
      majorant_bound_normal(maj, pf, Q0, αs, curτ, curN, n - 1, p)
      curτ += μ
      majorant_bound_special(maj, pf, Q0, μ, curτ, n, p)
      curN = n + 1
    end
  end

  majorant_bound_normal(maj, pf, Q0, αs, curτ, curN, fmpz(-1), p)

  # 1/Px[1+r] is majorized by                abs(1/c)
  #                             ----------------------------------
  #                             (1-x)^|k1-k2| * (1-x^2)^min(k1,k2)
  return (inv(abs(nacb(c, p), p), p), abs(k1-k2), min(k1,k2), maj, curτ)
end

function majorant_bound_normal(
  maj::hyp_majorant{narb},
  pf::hyp_majorant{S},
  Q0::S,
  αs::T,
  τ::Int,
  n0::fmpz, n1::fmpz,
  p::Int) where {S, T}

  for (i, j) in zip(coeffs(maj), coeffs(pf))
    max!(i, i, fraction_bound_normal(j, Q0, αs, τ, n0, n1, p), p)
  end
end

function majorant_bound_special(
  maj::hyp_majorant{narb},
  pf::hyp_majorant{S},
  Q0::S,
  μ::Int,
  τ::Int,
  n::fmpz,
  p::Int) where {S, T}

  for (i, j) in zip(coeffs(maj), coeffs(pf))
    max!(i, i, fraction_bound_special(j, Q0, μ, τ, n, p), p)
  end
end

# return upperbound on |a|
function abs_upper(a::nacb, p)
  mag = narb()
  t = mag_struct(0, 0)
  ccall((:mag_init, libarb), Nothing,
        (Ref{mag_struct},),
        t)
  ccall((:acb_get_mag, libarb), Nothing,
        (Ref{mag_struct}, Ref{nacb}),
        t, a)
  ccall((:arb_set_interval_mag, libarb), Nothing,
        (Ref{narb}, Ref{mag_struct}, Ref{mag_struct}, Int),
        mag, t, t, p)
  ccall((:mag_clear, libarb), Nothing,
        (Ref{mag_struct},),
        t)
  return mag
end


# u is a sxτxν array of the last s coefficients
# return sxτxν array of the "normalized residues"
function q_residual(Pθ, Px, u, λ::T, N::Int, τ::Int, ν::Int, Φ::Int) where T
  Fθ = parent(Pθ[1])
  F = base_ring(Fθ)
  θ = gen(Fθ)
  c = coeff(Px[end], 0)
  s = length(Pθ)-1
  @assert length(u) == s
  @assert length(u[1]) == τ
  @assert length(u[1][1]) == ν
  b = T[coeff(Pθ[1+i](λ+N+j+θ), k) for i in 0:s, j in 0:s, k in 0:τ-1]
  q = T[F() for i in 1:s, j in 1:τ, l in 1:ν]
  Q = nacb[nacb() for i in 1:s, j in 1:τ, l in 1:ν]
  v = T[F() for i in 1:s, j in 1:τ, l in 1:ν]
  for j in 0:s-1, k in τ-1:-1:0, l in 1:ν
    v[1+j,1+k,l] = zero(F)
    for jp in 1:s-j, kp in 0:τ-1-k
      v[1+j,1+k,l] += b[1+j+jp,1+j,1+kp]*u[jp][1+k+kp][l]
    end
    q[1+j,1+k,l] = c*v[1+j,1+k,l]
    for kp in 1:τ-1-k
      q[1+j,1+k,l] -= b[1+0,1+j,1+kp]*q[1+j,1+k+kp,l]
    end
    q[1+j,1+k,l] //= b[1+0,1+j,1+0]
    Q[1+j,1+k,l] = nacb(q[1+j,1+k,l], Φ)
  end
  return Q
end

# u is a sxτxν array of the last s coefficients
# return sxτxν array of the "normalized residues"
function q_residual(
  Pθ, Px,
  u::Array{RC, 3},
  λ::T,
  N::Int,
  τ::Int, ν::Int,
  Φ::Int) where {T, RC <: Union{narb, nacb}}

  Fθ = parent(Pθ[1])
  F = base_ring(Fθ)
  θ = gen(Fθ)
  s = length(Pθ)-1
  @assert (s, τ, ν) == size(u)
  c = RC(coeff(Px[end], 0), Φ)
  b = RC[RC(coeff(Pθ[1+i](λ+N+j+θ), k), Φ) for i in 0:s, j in 0:s, k in 0:τ-1]
  q = RC[RC() for i in 1:s, j in 1:τ, l in 1:ν]
  v = RC[RC() for i in 1:s, j in 1:τ, l in 1:ν]
  for j in 0:s-1, k in τ-1:-1:0, l in 1:ν
    zero!(v[1+j,1+k,l])
    for jp in 1:s-j, kp in 0:τ-1-k
      addmul!(v[1+j,1+k,l], b[1+j+jp,1+j,1+kp], u[jp,1+k+kp,l], Φ)
    end
    mul!(q[1+j,1+k,l], c, v[1+j,1+k,l], Φ)
    for kp in 1:τ-1-k
      submul!(q[1+j,1+k,l], b[1+0,1+j,1+kp], q[1+j,1+k+kp,l], Φ)
    end
    div!(q[1+j,1+k,l], q[1+j,1+k,l], b[1+0,1+j,1+0], Φ)
  end
  return q
end

function tail_bound(Pθ, Px, u, ns, αs, Z::zBranchInfo, δ, N, τ, ν, Φ)

  s = length(Pθ) - 1

  # for the residual it is required that none of λ+N, λ+N+1, ..., λ+N+s-1 are inidicial roots
  @assert isempty(ns) || N+s <= ns[1]

  # c/((1-z)^k1*(1-z^2)^k2) is supposed to majorize 1/Px[1+r], r = deg_θ(P)
  # exp(maj(z)) is the hhat(z)
  # finalτ is strict bound on the power of log(z) in the real solution f(z)
  (c, k1, k2, maj, finalτ) = equ_bound(Pθ, Px, Z.λ, ns, αs, τ, N, Φ)

  q = q_residual(Pθ, Px, u, Z.λ, N, τ, ν, Φ)

  Er = [zero(narb_poly) for l in 1:ν]

  if all(iszero, q)
    return Er
  end

  zmag = abs_upper(get_approx(Z, Φ), Φ)

  if (k1 > 0 || k2 > 0) && !isnegative(sub(zmag, 1, Φ))
    # don't know if the series converges
    for l in 1:ν, i in 0:δ-1
      setcoeff!(Er[l], i, pos_inf!(narb()))
    end
    return Er
  end

  for l in 1:ν
    f = zero(narb_poly)
    for i in 0:s-1
      # first take the max's of coeffs of log(z)^j/j!
      m = zero(narb)
      for j in 0:τ-1
        max!(m, m, abs(mul(q[1+i,1+j,l], N+i, Φ), Φ), Φ)
      end
      setcoeff!(f, i, m)
    end
    g = mullow(f, exp_series(neg!(series(maj, s, Φ)), s, Φ), s, Φ)
    for i in 0:s-1
      setcoeff!(g, i, max(zero(narb), div(coeff(g, i), N+i, Φ), Φ))
    end

    # f(z) is given by Σ_{i,j} u_{i,j}*z^(λ+i)*log(z)^j/j!
    # For each j, the remainder Σ_{i>=N} u_{i,j}*z^i is majorized by z^N*B(z)
    #     B(z) = g(z)*exp(maj(z))*c/((1-z)^k1*(1-z^2)^k2)
    #     TODO!!  when there are large initial roots so that ns is still
    #             not empty increase g(z) accordingly. for now we just
    @assert isempty(ns)

    # z^(λ+N-δ)
    zeval = narb_poly(zmag, 1)
    # TODO nonreal λ will require more schenanigans
    Er[l] = pow_series(zeval, narb(Z.λ + (N - δ), Φ), δ, Φ)
    for d in 0:δ-1
      setcoeff!(Er[l], d, abs!(coeff(Er[l], d)))
    end

    # c/((1-z)^k1*(1-z^2)^k2)
    mul!(Er[l], Er[l], c, Φ)
    f = pow_series(inv_series(sub(1, zeval, Φ), δ, Φ), narb(k1), δ, Φ)
    Er[l] = mullow(Er[l], f, δ, Φ)
    z2eval = mullow(zeval, zeval, δ, Φ)
    f = pow_series(inv_series(sub(1, z2eval, Φ), δ, Φ), narb(k2), δ, Φ)
    Er[l] = mullow(Er[l], f, δ, Φ)

    # g(z)
    t = one(narb_poly)
    f = narb_poly(coeff(g, 0))
    for i in 1:s-1
      t = mullow(t, zeval, δ, Φ)
      add!(f, f, mul(t, coeff(g, i), Φ), Φ)
    end
    Er[l] = mullow(Er[l], f, δ, Φ)

    # exp(maj(z))
    f = exp_series(eval_series(maj, zeval, δ, Φ), δ, Φ)
    Er[l] = mullow(Er[l], f, δ, Φ)

    # Σ_{j<finalτ} z^δ*log(z)^j/j!
    Er[l] = mullow(Er[l], logzpoly_series(Z, δ, finalτ, Φ), δ, Φ)
  end

  return Er
end

function get_multiplicity!(ns::Vector{fmpz}, iv::Vector{Vector{T}}, n::Int) where T
  @assert length(ns) == length(iv)
  nextu = Vector{T}[]
  while !isempty(ns) && ns[1] == n
    popfirst!(ns)
    push!(nextu, popfirst!(iv))
  end
  return nextu
end

function eval_log_poly(u::nacb_poly, τ::Int, logz::nacb, Φ::Int)
  i = τ-1
  t = coeff(u, τ-1-i)
  t2 = nacb();
  while (i-=1)≥0
    mul!(t, t, div!(t2, logz, i + 1, Φ), Φ)
    add!(t, t, coeff(u, τ-1-i), Φ)
  end
  return t
end


# evaluate the next coefficient u_N and update the arguments
# currently done in the field F and subject to blowup - do not use too much!
function next_u_coeffs!(u, τ, Pθ, λ, ns, iv, N, ν)
  s = length(Pθ)-1
  Fθ = parent(Pθ[1])
  F = base_ring(Fθ)
  θ = gen(Fθ)

  Pn = map(a -> a(θ + (λ + N)), Pθ)
  un = get_multiplicity!(ns, iv, N)
  rhs = [[zero(F) for k in 1:ν] for j in 1:τ]
  for i in 1:s
    # Pn[1+i] is a polynomial in θ. Apply it to u[i] where θ is
    # the shift operator and sub result from rhs
    for k in 0:degree(Pn[1+i])
      for j in 1:length(u[i])-k, l in 1:ν
        sub!(rhs[j][l], rhs[j][l], coeff(Pn[1+i], k)*u[i][j + k][l])
      end
    end
  end
  # μ = number of initial conditions at λ + N = multiplicity of root
  μ = length(un)
  for i in 0:μ-1
    @assert iszero(coeff(Pn[1], i))
  end
  for j in 1:τ
    push!(un, rhs[j])
  end
  for i in 1:τ
    for j in 1:i-1, l in 1:ν
      un[1+μ+τ-i][l] -= coeff(Pn[1], μ+j)*un[1+μ+τ-i+j][l]
    end
    for l in 1:ν
      un[1+μ+τ-i][l] //= coeff(Pn[1], μ)
    end
  end
  # trim zeros off un
  while length(un) > τ && all(iszero, un[end])
    pop!(un)
  end

  τ = max(τ, length(un))
  pop!(u)
  pushfirst!(u, un)
  return (un, τ)
end

function u_derivative!(uN::Vector{Vector{T}}, fac::T, τ::Int, ν::Int) where T
  for j in 0:τ-1, l in 1:ν
    mul!(uN[1+j][l], uN[1+j][l], fac)
    if j+1 < τ
      add!(uN[1+j][l], uN[1+j][l], uN[1+j+1][l])
    end
  end
end

function eval_basis_forward(
  Pθ::Vector{S}, Px::Vector{S},
  Z::zBranchInfo{T},
  ns::Vector{fmpz},
  αs::Vector{T},
  δ::Int,
  p::Int) where {S, T}

  p += 20

  ν = length(ns)    # number of initial conditions
  τ = 0             # strict bound on the power of log(z) thus far
  maxN = 10+20*p    # max number of terms to sum

  @assert length(ns) > 0
  @assert ns[1] == 0

  s = length(Pθ) - 1
  @assert s > 0

  Fθ = parent(Pθ[1])
  θ = gen(Fθ)
  F = base_ring(Fθ)

  # u is an array of the last s solutions
  # each solution is an array of coefficients of log^k(z)/k!  0 <= k < τ
  # each coefficient is an array of ν elements of F
  u = [Vector{T}[] for i in 1:s]

  # M is the matrix answer
  M = nacb_mat(δ, ν)

  # the ns will be consumed as we sum past them
  ns = deepcopy(ns)
  iv = [[i == j ? one(F) : zero(F) for i in 1:ν] for j in 1:ν]

  t = nacb()

  changedcount = 0
  N = 0

  while true
    may_stop = isempty(ns) #|| N+s < ns[1]
    if N > maxN
      error("internal error: parameters too large")
    end

    (uN, τ) = next_u_coeffs!(u, τ, Pθ, Z.λ, ns, iv, N, ν)

    # add uN*z^N to sum
    changed = false
    uN = deepcopy(uN)
    for d in 0:δ-1
      d == 0 || u_derivative!(uN, Z.λ+(N-d+1), τ, ν)
      for l in 1:ν
        zero!(t)
        for j in 0:τ-1
          if !iszero(uN[1+j][l])
            add!(t, t, mul(zpowlogpow(Z, N-d, j, p), uN[1+j][l], p), p)
          end
        end
        add!(t, t, M[1+d,l], p)
        changed = changed || !overlaps(M[1+d,l], t)
        M[1+d,l] = t
      end
    end

    N += 1

    if !may_stop
      continue
    end

    if !changed
      if (changedcount += 1) > 10
        break
      end
    else
      changedcount = 0
    end
  end

#println("forward used $N terms")

  Er = tail_bound(Pθ, Px, u, ns, αs, Z, δ, N, τ, ν, p)
  for d in 0:δ-1, l in 1:ν
    t = M[1+d,l]
    er = mul(coeff(Er[l], d), factorial(d, p), p)
    add_error!(t, er)
    M[1+d,l] = t
  end

  return M
end


function overlap_dist(a::Array{nacb, 3}, Δ::Array{nacb, 3}, b::Array{nacb, 3}, Φ::Int)
  (m, n, τ) = size(a)
  @assert (m, n, τ) == size(Δ) == size(b)
  dist = 0
  for i in 1:m, j in 1:n, k in 1:τ
    dist = max(dist, overlap_dist(a[i,j,k], Δ[i,j,k], b[i,j,k], Φ))
  end
  return dist
end

function overlap_dist(a::nacb, Δ::nacb, b::nacb, Φ::Int)
  if overlaps(a, b)
    return 0
  end
  q = div(magnitude(Δ), min(magnitude(a), magnitude(b)))
  if isspecial(q)
    return iszero(q) ? 0 : Φ
  else
    dist = exponent(q) + Φ
    if dist < 0
      return 0
    elseif dist > Φ
      return Φ
    else
      return Int(dist)
    end
  end
end


function eval_basis_bs(
  Pθ::Vector, Px::Vector,
  Z::zBranchInfo{T},
  ns::Vector{fmpz},
  αs::Vector{T},
  δ::Int,
  Φ::Int) where T

  Φ += 20

  ν = length(ns)    # number of initial conditions
  τ = 0             # strict bound on the power of log(z) thus far
  maxN = 100+20*Φ   # max number of terms to sum

  @assert length(ns) > 0
  @assert ns[1] == 0

  s = length(Pθ) - 1
  @assert s > 0

  Fθ = parent(Pθ[1])
  θ = gen(Fθ)
  F = base_ring(Fθ)
  @assert T == elem_type(F)

  # u is an array of the last s solutions
  # each solution is an array of coefficients of log^k(z)/k!  0 <= k < τ
  # each coefficient is an array of ν elements of F
  u = [Vector{T}[] for i in 1:s]

  # M is the matrix answer
  M = nacb_mat(δ, ν)

  # the ns will be consumed as we sum past them
  ns = deepcopy(ns)
  iv = [[i == j ? one(F) : zero(F) for i in 1:ν] for j in 1:ν]

  t = nacb()
  zn = one(nacb)

  N = 0
  N0 = -1
  uN0::Array{nacb, 3} = [zero(nacb) for i in 1:s, j in 1:ν, k in 1:1]
  M0 = nacb_mat(δ, ν)

  while true
    if N > maxN
      error("internal error: parameters too large")
    end

    (uN, τ) = next_u_coeffs!(u, τ, Pθ, Z.λ, ns, iv, N, ν)

    # add un*z^N to sum
    for d in 0:δ-1
      if d > 0
        if d == 1
          uN = deepcopy(uN)
        end
        u_derivative!(uN, Z.λ+(N-d+1), τ, ν)
      end
      for l in 1:ν
        zero!(t)
        for j in 0:τ-1
          if !iszero(uN[1+j][l])
            add!(t, t, mul(zpowlogpow(Z, N-d, j, Φ), uN[1+j][l], Φ), Φ)
          end
        end
        add!(t, t, M[1+d,l], Φ)
        setindex!(M, t, 1+d,l)
      end
    end

    N += 1

    if isempty(ns) && N≥δ
      # for u_{N0-1},...,u_{N0-s}, map the 3D array u of F elem's
      # into 2D array uN0 of poly's in Λ (it‘s still 3D)
      N0 = N
      if τ != 1
        uN0 = [zero(nacb) for i in 1:s, j in 1:ν, k in 1:τ]
      end
      for i in 1:s, l in 1:ν
        for j in 0:τ-1
          if j < length(u[i])
            uN0[i,l,1+τ-1-j] = nacb(u[i][1+j][l], Φ)
          else
            uN0[i,l,1+τ-1-j] = zero(nacb)
          end
        end
      end
      break
    end
  end

  # set up sum for [N0, N1)
  numeq = elem_type(Fθ)[zero(Fθ) for j in 1:s, i in 0:τ-1]
  Pn = map(a -> a(θ + Z.λ), Pθ)
  den = Pn[1]^(τ-1)
  for j in 1:s
    d = -Pn[1+j]
    g = den
    numeq[j,1+0] = g*d
    for i in 1:τ-1
      d = divexact(derivative(d)*Pn[1], i) - d*derivative(Pn[1])
      g = divexact(g, Pn[1])
      numeq[j,1+i] = d*g
    end
  end
  den *= Pn[1]

  # sum [N0, N1)
  @assert τ>0
  @assert N0≥δ
  N1 = N0+10
  S = sum_bs_start(numeq, den, Z.biz, Z.λ, δ, N1, N0, Φ)
  overlap_count = 0
  d1 = Φ
  NΔ = N1
  while N1 < maxN
    N2 = N1 + clamp(NΔ, fld(N1+100,64), N1)
    d2 = sum_bs_continue(S, N2, Φ)
    if d2 < 1
      overlap_count += 1
      NΔ = 10
    else
      overlap_count = 0
      NΔ = d1-d2 <= 0 ? N1 : trunc(Int, (N2-N1)*0.25*(d2/(d1-d2)))
    end
    N1 = N2
    overlap_count < 2 || break
  end

  # accumulate sums for [N0, N1) and add to sum for [0,N0)
  Σ1a = get_sum(S, Φ)
  fdz = [zero(nacb) for k in 1:τ]
  for d in 0:δ-1, l in 1:ν
    # fdz = z^(d-N0-λ)*f^(d)(z) for the l^th initial value
    for i in 1:s, k in 0:τ-1, h in 0:k
      if i == 1 && h == 0
        mul!(fdz[1+k], Σ1a[1+d,i,1+h], uN0[i,l,1+k-h], Φ)
      else
        mul!(t, Σ1a[1+d,i,1+h], uN0[i,l,1+k-h], Φ)
        add!(fdz[1+k], fdz[1+k], t, Φ)
      end
    end
    M[1+d,l] = add(M[1+d,l], zpowlogpoly_eval(Z, N0-d, fdz, Φ), Φ)
  end

  # compute u_{N1-1},...,u_{N1-s}
  uN1 = mullow(get_approx(S.M, Φ), uN0, Φ)
  uN1t = nacb[uN1[i,l,1+τ-j] for i in 1:s, j in 1:τ, l in 1:ν]

  # bound sum for [N1, ∞) and add to sum for [0,N1)
  Er = tail_bound(Pθ, Px, uN1t, ns, αs, Z, δ, N1, τ, ν, Φ)
  for d in 0:δ-1, l in 1:ν
    getindex!(t, M, 1+d,l)
    add_error!(t, mul(coeff(Er[l], d), factorial(ZZ(d)), Φ))
    setindex!(M, t, 1+d,l)
  end

  return M
end


# Return δ by ν matrix M such that for any vector iv of ν initial values the
# solution f(z) determined by these initial values and its derivatives is M.iv.
#  f(z) = Σ_{i,j} z^(λ-0+i)*log(z)^j/j!*( u[i,j] )
# f'(z) = Σ_{i,j} z^(λ-1+i)*log(z)^j/j!*( (λ+i)*u[i,j] + u[i,j+1] )
function eval_basis(
  Pθ::Vector, Px::Vector, # equation with θ on the inside and out, but always on the left
  Z::zBranchInfo,         # λ, z and log(z)
  ns::Vector{fmpz},       # indicial roots λ+ns[1], λ+ns[2], ..., λ+ns[ν]
  αs::Vector{T},          # all indicial roots
  δ::Int,                 # calculate f(z), f'(z), ... f^(δ-1)(z)
  Φ::Int) where T         # precision

  M1 = eval_basis_bs(Pθ, Px, Z, ns, αs, δ, Φ)
  if @debug()
    M2 = eval_basis_forward(Pθ, Px, Z, ns, αs, δ, Φ)
    for i in 1:δ, j in 1:length(ns)
      @assert overlaps(M1[i,j], M2[i,j])
    end
  end
  return M1
end


function eval_bases(Pθ, Px, z, ρ, αs::Vector{T}, δ, Φ) where T
  M = eval_basis(Pθ, Px, zBranchInfo{T}(ρ[1][1], z), ρ[1][2], αs, δ, Φ)
  for i in 2:length(ρ)
    M = hcat(M, eval_basis(Pθ, Px, zBranchInfo{T}(ρ[i][1], z), ρ[i][2], αs, δ, Φ))
  end
  return M
end


#### nFn-1 evaluation with parameters in a field F #############################

# return array of (λ_i, [n_i1 = 0, n_i2, ...], Ind_i) such that
#   the λ_i are distinct mod 1
#   0 = n_i1 <= n_i2 <= n_i3 .... <= n_ik are integers
#   for each i, the unordered list {λ_i + n_ij}_j = the unordered list {a[Ii[j]]}_j
#   the Ind_i[j] are of course all distinct
function partition_mod1(a::Vector{T}) where T
  F = parent(a[1])
  r = Tuple{elem_type(F), Vector{elem_type(ZZ)}, Vector{Int}}[]
  for i in 1:length(a)
    found = false
    for j in 1:length(r)
      ok, k = isinteger_with_integer(a[i] - r[j][1])
      if ok
        push!(r[j][2], k)
        push!(r[j][3], i)
        found = true
        break
      end
    end
    if !found
      push!(r, (a[i], [ZZ(0)], [i]))
    end
  end
  for j in 1:length(r)
    sort!(r[j][2])
    r[j] = (r[j][1] + r[j][2][1], r[j][2] .- r[j][2][1], r[j][3]) 
  end
  return r
end

# series in a + b*ε
function rgamma_series(a::nacb, b::Int, ord::Int, p::Int)
  z = nacb_poly(a, b)
  ccall((:acb_poly_rgamma_series, libarb), Nothing,
        (Ref{nacb_poly}, Ref{nacb_poly}, Int, Int),
        z, z, ord, p)
  return z
end

function gamma_series(a::nacb, b::Int, ord::Int, p::Int)
  z = nacb_poly(a, b)
  ccall((:acb_poly_gamma_series, libarb), Nothing,
        (Ref{nacb_poly}, Ref{nacb_poly}, Int, Int),
        z, z, ord, p)
  return z
end

function rising_factorial_series(a::nacb, b::Int, t::fmpz, ord::Int, p::Int)
  @assert t >= 0
  z = nacb_poly(a, b)
  # TODO handle large t
  ccall((:acb_poly_rising_ui_series, libarb), Nothing,
        (Ref{nacb_poly}, Ref{nacb_poly}, UInt, Int, Int),
        z, z, UInt(t), ord, p)
  return z
end

function pi!(z::narb, p::Int)
  ccall((:arb_const_pi, libarb), Nothing,
        (Ref{narb}, Int),
        z, p)
  return z
end

function pi!(z::nacb, p::Int)
  GC.@preserve z begin
    t = ccall((:acb_real_ptr, libarb), Ptr{narb}, (Ref{nacb},), z)
    ccall((:arb_const_pi, libarb), Nothing, (Ptr{narb}, Int), t, p)
    t = ccall((:acb_imag_ptr, libarb), Ptr{narb}, (Ref{nacb},), z)
    ccall((:arb_zero, libarb), Nothing, (Ptr{narb},), t)
  end
  return z
end

# π*ε*csc(x*π*ε)
function cscpi_series(x::fmpz, ord::Int, p::Int)
  z = nacb_poly()
  setcoeff!(z, 0, nacb(1//x, p))
  t = narb()
  for n in 1:div(ord-1,2)
    pi!(t, p)
    pow!(t, t, UInt(2*n), p)
    mul!(t, t, divexact(bernoulli(2*n),factorial(ZZ(2*n)))*
               (-1)^(n+1)*(ZZ(4)^n-2)*x^(2*n-1), p) # TODO
    setcoeff!(z, 2*n, nacb(t))
  end
  return z
end


# integers am[1] ≤ am[2] ≤ ... ≤ am[m]
# the numerator parameters are λ + am[1], ..., λ + am[m], as[1], ... , as[.]
# produce the initial values for the exponents λ + am[1], ..., λ + am[m]
# computed by replacing λ + am[1+i] -> λ + am[1+i] + i*ε and letting ε -> 0
# no difference λ - as[i] is an integer
#=
  keep in mind:
    g(a,ε,t) := Γ(a+ε)/(1-(a+ε))_t
              = (-1)^t*Γ(a-t+ε)               for integer t, nice for a > t
              = (-1)^a*π*csc(π*ε)/Γ(1-a+t-ε)  for integer a, nice for a ≤ t

  hat(a)[k] is the a-vector with the k^th entry omitted

  derivation:  of concern in the formula

      a            Γ(b)*Γ(hat(a)[k] - a[k])  a[k]   a[k], 1+a[k]-b       p-q
  pFq( ;-1/z) =  Σ ------------------------ z    F(                ; (-1)   *z)
      b         k≤p Γ(b-a[k])*Γ(hat(a)[k])         1+a[k]-hat(a)[k]

  is the coefficient of z^(λ + n)*ε^0 as a polynomial in log(z). We don't care
  about negative powers of ε because the total contribution of them is zero,
  and we don't care about positive powers of ε because ε -> 0.

  recall that the a vector is partition into
    {λ + am[1+i]}_{0≤i<m}, and the rest as[1], ... , as[p-m]

  the coefficient of z^(λ + n) is that of z^(λ + n) in

        Γ(b) * (Π_{i<m,i≠k} Γ(a[1+i]-a[1+k])) * Γ(as-a[1+k])
    Σ   --------------------------------------------------------
  0≤k<m     Γ(b-a[1+k]) * (Π_{i<m,i≠k} Γ(a[1+i])) * Γ(as)

      * exp(a[1+k]*log(z))

                       (a[1+k])_t * (1+a[1+k]-b)_t                  (-1)^(t*(p-q))*z^t
      * Σ  ------------------------------------------------------- -------------
       0≤t (Π_{i<m,i≠k} (1+a[1+k]-a[1+i])_t) * (1+a[1+k]-as)_t       t!

  using a[1+i] = λ+am[1+i]+i*ε for i<m

    Γ(b)             Π_{i<m,i≠k} Γ(am[1+i]-am[1+k]+(i-k)*ε)
  = ---- *   Σ    Σ  ----------------------------------------------
    Γ(as)  0≤k<m 0≤t Π_{i<m,i≠k} (1+am[1+k]-am[1+i]+(k-i)*ε)_t

                (λ+am[1+k]+k*ε)_t              * (1+λ+am[1+k]+k*ε-b)_t
          * ------------------------------------------------------------------
            (Π_{i<m,i≠k} Γ(λ+am[1+i]+i*ε)) * Γ(b-λ-am[1+k]-k*ε)

            Γ(as-λ-am[1+k]-k*ε)                    (-1)^(t*(p-q))*z^(t+λ+am[1+k])
          * --------------------- exp(k*ε*log(z)) -------------------------------
            (1+λ+am[1+k]+k*ε-as)_t                            t!

 ...

    Γ(b)                                    g(am[1+i]-am[1+k], (i-k)*ε, t)/(-1)^t
  = ---- *   Σ    Σ (λ+am[1+k]+k*ε)_t    Π  -------------------------------------
    Γ(as)  0≤k<m 0≤t                  i<m,i≠k           Γ(λ+am[1+i]+i*ε)

                    Γ(as-λ-am[1+k]-t-k*ε)
                  * ---------------------
                    Γ(b-λ-am[1+k]-t-k*ε)

                       (k*ε)^i*log(z)^i z^(λ+t+am[1+k])
                  * Σ  ---------------- ---------------
                   0≤i         i!             t!
=#
function hyp_initial_values_at_∞(
  λ::T, am::Vector{fmpz}, as::Vector{T},
  b::Vector{T},
  p::Int) where T

  p += length(am) + length(as)

  m = length(am)
  r = [zero(nacb) for j in 1:m]
  j = 0
  while j < m
    n = am[1+j]
    # μ = multiplicity of exponent λ + n
    μ = 1
    while j+μ < m && am[1+j+μ] == n
      μ += 1
    end
    # compute coeffs r[1+j], ..., r[1+j+μ-1]
    # of z^(λ+n), ..., z^(λ+n)*log^(μ-1)(z)/(μ-1)!
    for k in 0:m-1
      # t is the index of the term in the (q+1)F(p-1) series
      t = n - am[1+k]
      if t < 0
        continue
      end
      # compute the order required for the power series O(ε^ord)
      ord = 1
      for i in 0:m-1
        if i != k && am[1+i]-am[1+k]-t <= 0
          # a gamma factor has a simple pole
          ord += 1
        end
      end
      s = one(nacb_poly)
      # handle integer differences between numerator parameters
      for i in 0:m-1
        if i == k
          s2 = rising_factorial_series(nacb(λ+am[1+i], p), i, t, ord, p)
          s = mullow(s, s2, ord, p)
        else
          s2 = rgamma_series(nacb(λ+am[1+i], p), i, ord, p)
          s = mullow(s, s2, ord, p)
          d = am[1+i]-am[1+k]-t
          if d > 0
            s2 = gamma_series(nacb(d, p), i-k, ord, p)
            s = mullow(s, s2, ord, p)
          else
            s2 = rgamma_series(nacb(1-d, p), k-i, ord, p)
            s2 = mullow(s2, cscpi_series(ZZ(i-k), ord, p), ord, p)
            s = mullow(s, s2, ord, p)
            isodd(d) && neg!(s)
          end
        end
      end
      # no difference A - λ - am[1 + i] is an integer
      for A in as
        s2 = gamma_series(nacb(A-λ-am[1+k]-t, p), -k, ord, p)
        s = mullow(s, s2, ord, p)
      end
      # the difference B - λ - am[1 + i] could be an integer, but it's rgamma
      for B in b
        s2 = rgamma_series(nacb(B-λ-am[1+k]-t, p), -k, ord, p)
        s = mullow(s, s2, ord, p)
      end
      # now determine r[1+j], ..., r[1+j+μ-1]. They are already zeroed.
      f = (iseven(t) ? 1 : -1)//factorial(t)
      for i in 0:min(ord, μ)-1
        add!(r[1+j+i], r[1+j+i], mul(coeff(s, ord-1-i), f, p), p)
        mul!(f, f, k)
      end
    end
    j += μ
  end
  f = one(nacb)
  for A in as
    mul!(f, f, rgamma(nacb(A, p), p), p)
  end
  for B in b
    mul!(f, f, gamma(nacb(B, p), p), p)
  end
  iv = nacb_mat(m, 1)
  for i in 1:m
    iv[i,1] = mul(r[i], f, p)
  end
  return iv
end

function compute_f_near_0(a::Vector{T}, b::Vector{T}, z::nacb, Φ::Int) where T
  F = parent(a[1])
  Fx, x = PolynomialRing(F, "x")
  # 0 -> z
  Pθ, Px = hyp_equ(a, b, Fx(0)//Fx(1), x//Fx(1))
  αs = push!(F(1) .- b, F(0))
  e = eval_basis(Pθ, Px, zBranchInfo{T}(F(0), z), [ZZ(0)], αs, 1, Φ)
  return e[1,1]
end

function compute_f_near_∞(a::Vector{T}, b::Vector{T}, z::nacb, p::Int) where T
  F = parent(a[1])
  Fx, x = PolynomialRing(F, "x")
  # ∞ -> z
  Pθ, Px = hyp_equ(a, b, Fx(0)//Fx(1), -Fx(1)//x)
  αs = a
  nz = neg(z)
  rz = inv(nz, p)
  s = zero(nacb)
  for (λ, ns, indices) in partition_mod1(a)
    arest = T[a[j] for j in 1:length(a) if !(j in indices)]
    iv = hyp_initial_values_at_∞(λ, ns, arest, b, p)
    e = eval_basis(Pθ, Px, zBranchInfo{T}(λ, rz, nz), ns, αs, 1, p)
    add!(s, s, mul(e, iv, p)[1,1], p)
  end
  return s
end

# recurrence eq for buehring's coeffs Ak(a+s,b+s)
function buehring_equ(a::Vector{T}, b::Vector{T}, s::Int, θ) where T
  Fθ = parent(θ)
  t::T = 1-s-a[1]
  eq = [θ, -(θ-1+t)*(θ-1+b[1]-a[1])]
  for i in 2:length(b)
    t += b[i-1]-a[i]
    eq = hadamard_product(eq, [(θ-1+t), -one(Fθ)])
    # before cauchy uses the diff equ, ensure that eq matches the generic case
    if degree(eq[1]) < i
      eq = eq .* divexact(θ*prod(θ+t-j for j in 1:i-1), eq[1])
    end
    eq = cauchy_product(eq, [θ, -(θ-1+b[i]-a[i])])
    eq = hadamard_product(eq, [one(Fθ), -(θ-1+t)])
  end
  return eq
end

function const_euler_gamma(p::Int)
  z = narb()
  ccall((:arb_const_euler, libarb), Nothing,
        (Ref{narb}, Int),
        z, p)
  return z
end

function rgamma_jet2(a::T, p::Int) where T
  z = rgamma_series(nacb(a, p), 1, 2, p)
  return (coeff(z,0), coeff(z,1))
end

function rising_factorial_jet2(a::T, n::fmpz, p::Int) where T
  z = rising_factorial_series(nacb(a, p), 1, n, 2, p)
  return (coeff(z,0), coeff(z,1))
end


# for intcase
function _hyp2f1_initial_values_at_1!(iv::nacb_mat, a::T, b::T, c::T, σ::fmpz, Φ::Int) where T
  r1 = const_euler_gamma(Φ)
  rσa0, rσa1 = rgamma_jet2(σ+a, Φ)
  rσb0, rσb1 = rgamma_jet2(σ+b, Φ)
  c0 = mul(rσa0, rσb0, Φ)
  c1 = add(mul(rσa1, rσb0, Φ), mul(rσa0, rσb1, Φ), Φ)
  if iszero(σ)
    iv2 = neg(c0)
    iv1 = sub(c1, mul(ldexp!(c0, c0, 1), r1, Φ), Φ)
  else
    rab = mul(rgamma(a, Φ), rgamma(b, Φ), Φ)
    if σ > 0
      iv1 = mul(c0, factorial(σ-1), Φ)
      rσ0, rσ1 = rgamma_jet2(1+σ, Φ)
      fσa = rising_factorial(a, σ, Φ)
      fσb = rising_factorial(b, σ, Φ)
      t = sub(mul(c0, r1, Φ), c1, Φ)
      iv2 = add(mul(rσ1, rab, Φ), mul(mul(rσ0, t, Φ), mul(fσa, fσb, Φ), Φ), Φ)
    else; @assert σ < 0
      iv1 = mul(rab, factorial(-1-σ), Φ)
      rσ0, rσ1 = rgamma_jet2(1-σ, Φ)
      fσa0, fσa1 = rising_factorial_jet2(σ+a, -σ, Φ)
      fσb0, fσb1 = rising_factorial_jet2(σ+b, -σ, Φ)
      t = add(add(mul(fσa1, fσb0, Φ), mul(fσa0, fσb1, Φ), Φ), mul(mul(fσa0, fσb0, Φ), r1, Φ), Φ)
      t = mul(t, mul(rσ0, rab, Φ), Φ)
      iv2 = sub(add(t, mul(rσ1, c0, Φ), Φ), mul(rσ0, c1, Φ), Φ)
    end
    if iseven(σ)
      neg!(iv2)
    end
  end
  Γc = gamma(c, Φ)
  iv[1,1] = mul(iv1, Γc, Φ)
  iv[2,1] = mul(iv2, Γc, Φ)
  return iv
end

# in the case of non-integral σ, output its iv first
function hyp_initial_values_at_1(
  a::Vector{T}, b::Vector{T},
  Φ::Int,
  succ_method::Bool) where T

  q = length(b)
  @assert length(a) == q+1
  σ::T = sum(b) - sum(a)
  Φ += 20
  m = trunc(Int, Φ+q)
  N = trunc(Int, Φ+q)
  iv = nacb_mat(q+1, 1)

  if q < 1
    iv[1,1] = one(nacb)
    return iv
  end

  intcase, ZZσ = isinteger_with_integer(σ)
  if intcase
    @assert q == 1 # q > 1 not implemented too messy
    _hyp2f1_initial_values_at_1!(iv, a[1], a[2], b[1], ZZσ, Φ)
    return iv
  end

  Γb1 = gamma(b[1], Φ)
  Γmσ = gamma(-σ, Φ)

  Γboa = div(Γb1, gamma(a[1], Φ), Φ)
  for i in 2:q
    mul!(Γboa, Γboa, div(gamma(b[i], Φ), gamma(a[i], Φ), Φ), Φ)
  end

  # u0
  iv[1,1] = mul(div(Γmσ, gamma(a[end], Φ), Φ), Γboa, Φ)
  if q == 1
    t = div(Γb1, Γmσ, Φ)
    t = mul(t, div(neg!(pi!(narb(), Φ)), mul(sin_pi(σ, Φ), σ, Φ), Φ), Φ)
    iv[1,2] = mul(t, mul(rgamma(b[1]-a[1], Φ), rgamma(b[1]-a[2], Φ), Φ), Φ)
    return iv
  end

  Γfac = mul(Γboa, div(rising_factorial(a[end], m, Φ), gamma(σ+a[end]+m, Φ), Φ), Φ)

  F = parent(a[end])
  Fθ, θ = PolynomialRing(F, "θ")

  # the sum of the first m terms of the defining series
  Φp = Φ+4*nbits(m)
  S = sum_bs_start(fill(prod(θ-1+ai for ai in a),1,1), θ*prod(θ-1+bi for bi in b),
                   BiElem{T}(one(F), F), zero(F), q, m, 1, Φp)
  plain = get_sum(S, Φp)
  add!(plain[1,1,1], plain[1,1,1], 1, Φp)

  # the vi coeffs
  if succ_method
    for i in 0:q-1
      Aeq = buehring_equ(a, b, i, θ)
      eq = hadamard_product(Aeq, [(θ+σ-i)*(θ-1+σ+a[end]+m-i), -(θ-1+σ-i)])
      tail = sum_bs(N, eq, initial_values(eq, [F(1)]), Φ+length(eq)+4*nbits(N))
      i==0 || mul!(Γfac, Γfac, σ-i+a[end]+m, Φ)
      t2 = mul(mul(Γfac, inv(σ-i), Φ), tail, Φ)
      t1 = plain[1+i,1,1]
      iv[2+i,1] = div(add(t1, t2, Φ), (-1)^i*factorial(ZZ(i)), Φ)
    end
  else
    Aeq = buehring_equ(a, b, 0, θ)
    S = [zero(narb) for i in 0:q-1]
    for j in 0:q-1
      eq = hadamard_product(Aeq, [(θ+σ)*(θ-1+σ+a[end]+m), -(θ-1-j+σ)])
      S[1+j] = sum_bs(N, eq, initial_values(eq, [F(1)]), Φ+length(eq)+4*nbits(N))
    end
    for i in 0:q-1
      tf = inv(σ)
      tail = zero(narb)
      for j in 0:i
        add!(tail, tail, mul(S[1+j], (-1)^(i-j)*binomial(fmpz(m),fmpz(i-j))*tf, Φ), Φ)
        tf *= (a[end]+m+j)//(1+j-σ)
      end
      t2 = mul(tail, Γfac, Φ)
      t1 = div(plain[1+i,1,1], (-1)^i*factorial(ZZ(i)), Φ)
      iv[2+i,1] = add(t1, t2, Φ)
    end
  end
  return iv
end


function compute_f_near_1(a::Vector{T}, b::Vector{T}, z::nacb, p::Int) where T
  p += 20
  n = length(a)
  @assert n-1 == length(b)
  F = parent(a[1])
  Fx, x = PolynomialRing(F, "x")
  w = F(3//4)
  # 0 -> w
  Pθ, Px = hyp_equ(a, b, Fx(0)//Fx(1), x//Fx(1))
  αs = push!(F(1) .- b, F(0))
  e0 = eval_basis(Pθ, Px, zBranchInfo{T}(F(0), w), [ZZ(0)], αs, n, p)
  for i in 2:2:n
    e0[i,1] = neg!(e0[i,1]) # correct odd derivatives
  end
  # roots of indicial equation are 0:n-2 and σ
  σ::T = sum(b) - sum(a)
  αs = push!(F.(0:n-2), σ)
  intcase, ZZσ = isinteger_with_integer(σ)
  if n < 2
    ρ = [(σ, [ZZ(0)])]
  elseif intcase
    roots = sort!(push!(ZZ.(0:n-2), ZZσ))
    ρ = [(F(roots[1]), roots .- roots[1])]
  else
    ρ = [(σ, [ZZ(0)]), (F(0), ZZ.(0:n-2))]
  end
  # 1 -> w and 1 -> z
  Pθ, Px = hyp_equ(a, b, Fx(0)//Fx(1), (1-x)//Fx(1))
  e2 = eval_bases(Pθ, Px, 1-w, ρ, αs, n, p)
  # 0 -> w -> 1 -> z
  sol = solve(e2, e0, p)
  if @debug() && (length(b) < 2 || !intcase)
    sol2 = hyp_initial_values_at_1(a, b, p, false)
    for i in 1:length(a)
      @assert overlaps(sol2[i,1], sol[i,1])
    end
    sol3 = hyp_initial_values_at_1(a, b, p, true)
    for i in 1:length(a)
      @assert overlaps(sol3[i,1], sol[i,1])
    end
  end
  e1 = eval_bases(Pθ, Px, sub(1, z, p+60), ρ, αs, 1, p)
  return mul(e1, sol, p)[1,1]
end

function compute_f_anywhere(a::Vector{T}, b::Vector{T}, z::nacb, p::Int) where T
  F = parent(a[1])
  Fx, x = PolynomialRing(F, "x")
  # 0 -> (1-sqrt(1-z))/(1+sqrt(1-z))
  Pθ, Px = hyp_equ(a, b, -2*Fx(a[1])//(1+x), 4*x//(1+x)^2)
  αs = push!(F(1) .- b, F(0))
  t = sqrt(sub(1, z, p+60), p+60)
  s = add(1, t, p+60)
  w = div(sub(1, t, p+60), s, p)
  e = eval_basis(Pθ, Px, zBranchInfo{T}(F(0), w), [ZZ(0)], αs, 1, p)
  ldexp!(s, s, -1)
  return mul(pow(s, -2*a[1], p), e[1,1], p)
end

# The non-regularized function is undefined if any bi is in {0, -1, -2, ...}.
# The regularized function is not implemented and may be defined via the
# initial values at the exponents 0 and the 1 - bi:
#   - any non-integer exponent has value 0
#   - only the greatest integer exponent n has nonzero value and this value is
#      prod_i (ai)_n / prod_i Γ(bi + n) * z^n / n!  (with no log(z))
function compute_pfq(a::Vector{T}, b::Vector{T}, z::nacb, p::Int) where T
  Nemo.arb_check_precision(p+200)
  p += 20
  zz = convert(Complex{Float64}, z)
  if length(a) < length(b) + 1
    return compute_f_near_0(a, b, z, p)
  elseif length(a) > length(b) + 1
    contains_zero(z) && error("oops not implemented")
    return compute_f_near_∞(a, b, z, p)
  else
    if abs2(zz) < (7/8)^2
      return compute_f_near_0(a, b, z, p)
    elseif abs2(zz) > (8/7)^2
      return compute_f_near_∞(a, b, z, p)
    elseif abs2(zz) < 2*real(zz) - 0.75
      return compute_f_near_1(a, b, z, p)
    else
      return compute_f_anywhere(a, b, z, p)
    end
  end
end

#### tests #####################################################################

function θ3(a::acb, b::acb)
  R = parent(a)
  z = [R(), R(), R(), R()]
  ccall((:acb_modular_theta, libarb), Nothing,
        (Ref{acb}, Ref{acb}, Ref{acb}, Ref{acb}, Ref{acb}, Ref{acb}, Int),
        z[1], z[2], z[3], z[4], a, b, precision(R))
  return z[3]
end

function gamma(a::Rational, C::AcbField)
  R = ArbField(precision(C))
  return C(gamma(QQ(a), R))
end

function hypergeometric_pfq(a::Vector, b::Vector, z::acb)
  C = parent(z)
  A = fmpq[QQ(i) for i in a]
  B = fmpq[QQ(i) for i in b]
  return C(compute_pfq(A, B, nacb(z), precision(C)))
end


function run_tests()
println()
println("***************************** tests *********************************")

if @debug()
  @testset "initial values" begin
    CC = AcbField(100)
    for z in [CC(99//100), CC(101//100)], c in 1:5, k in 1:5
      @test isfinite(hypergeometric_pfq([1//3,2//3],[c],z))
      @test isfinite(hypergeometric_pfq([k+1//3,2//3],[c],z))
      @test isfinite(hypergeometric_pfq([1//3,k+2//3],[c],z))
    end
  end
end

@testset "singular points" begin
  CC = AcbField(53)

  z = CC("[0.0 +/- 1.0e-10]")
  @test overlaps(hypergeometric_pfq([1,2,3],[4,5],z), one(CC))

  z = zero(CC)
  @test overlaps(hypergeometric_pfq([1,2,3],[4,5],z), one(CC))

  z = CC("[1.0 +/- 1.0e-10]")
  @test overlaps(hypergeometric_pfq([1,2,3],[4,5],z), 120-12*const_pi(CC)^2)

  z = one(CC)
  @test overlaps(hypergeometric_pfq([1,2,3],[4,5],z), 120-12*const_pi(CC)^2)

  z = CC("[1.0 +/- 1.0e-10]")
  @test !isfinite(hypergeometric_pfq([5,4,3],[2,1],z))

  z = one(CC)
  @test !isfinite(hypergeometric_pfq([5,4,3],[2,1],z))
end

@testset "algebraic cases" begin
  CC = AcbField(100)

  for j in [CC(1//7+2000), CC(1//7-2000),
            CC(1//7+1728), CC(1728), CC(1//7-1728),
            CC(1//7+200), CC(1//7-200),
            CC(1//7+64), CC(-64), CC(1//7-64),
            CC(1//7+2), CC(1//7-2)]
    f1 = hypergeometric_pfq([-1//12,1//4],[2//3],1728/j)
    f2 = hypergeometric_pfq([1//4, 7//12],[4//3],1728/j)
    t = j/3^3*f1^3/f2^3
    @test overlaps(j, 27*t*(t+8)^3/(t-1)^3)

    f1 = hypergeometric_pfq([-1//24,7//24],[3//4],1728/j)
    f2 = hypergeometric_pfq([5//24,13//24],[5//4],1728/j)
    t = j/2^4*f1^4/f2^4
    @test overlaps(j, 16*(t^2+14*t+1)^3/(t*(t-1)^4))

    f1 = hypergeometric_pfq([-1//60,19//60],[4//5],1728/j)
    f2 = hypergeometric_pfq([11//60,31//60],[6//5],1728/j)
    t = j*f1^5/f2^5
    @test overlaps(j, (t^4+228*t^3+494*t^2-228*t+1)^3/(t*(t^2-11*t-1)^5))

    f1 = hypergeometric_pfq([-1//42,13//42,9//14],[4//7,6//7],1728/j)
    f2 = hypergeometric_pfq([5//42,19//42,11//14],[5//7,8//7],1728/j)
    f3 = hypergeometric_pfq([17//42,31//42,15//14],[9//7,10//7],1728/j)
    t = j*f1*f2^2/f3^3
    @test overlaps(j, (t^2-t+1)^3*(t^6+229*t^5+270*t^4-1695*t^3+1430*t^2-
                                      235*t+1)^3/((t^2-t)*(t^3-8*t^2+5*t+1)^7))

    f1 = hypergeometric_pfq([-1//10,1//10,1//2],[1//5,4//5],-64/j)
    f2 = hypergeometric_pfq([1//10,3//10,7//10],[2//5,6//5],-64/j)
    f3 = hypergeometric_pfq([7//10,9//10,13//10],[8//5,9//5],-64/j)
    t = j*f2^2/f3*(j*f1^2 + f2*f3)/(j*f2^3 + f1*f3^2)
    @test overlaps(j, (t^2-1)*(t^2-4*t-1)^5/(t*(t^2+t-1)^5))
  end
end

@testset "function zeros" begin
  CC = AcbField(100)
  i = onei(CC)
  ipi = const_pi(CC)*i

  for m in [CC(1//10), CC(2//3), CC(9//10)]
    f1 = hypergeometric_pfq([1//3,2//3],[1],m)
    f2 = hypergeometric_pfq([1//3,2//3],[1],1-m)
    f3 = hypergeometric_pfq([5//6,6//6,7//6],[3//2,3//2],1-m)
    τ = sqrt(CC(-1//3))*f2/f1
    z = (2 + 3*τ)/4 - sqrt(1-m)*f3/(6*ipi*f1)
    @test contains_zero(θ3(z,2*τ)*θ3(z,6*τ) + exp(2*ipi*τ)*θ3(z+τ,2*τ)*θ3(z-3*τ,6*τ))
  end

  for j in [CC(10000//3), CC(10//3)]
    f1 = hypergeometric_pfq([1//12,5//12],[1],1728/j)
    f2 = hypergeometric_pfq([1//12,5//12],[1//2],1-1728/j)
    f3 = hypergeometric_pfq([1//3,2//3,1],[3//4,5//4],1-1728/j)
    τ = i*(2*gamma(1//2,CC)/gamma(7//12,CC)/gamma(11//12,CC)*f2/f1 - 1)
    z = (1 + τ)/2 + sqrt(CC(2//3))/ipi*(1-1728/j)^(1//4)*f3/f1
    @test contains_zero(ellipwp(z, τ))
  end
end

@testset "entire case" begin
  CC = AcbField(100)

  for z in [CC(1), CC(-1), CC(10), CC(-10)]
    @test overlaps(hypergeometric_pfq([1],[2],z), (exp(z)-1)/z)
  end
end

@testset "divergent case" begin
  CC = AcbField(100)

  # 0! - 1! + 2! - 3! + ... = 0.59...
  @test overlaps(hypergeometric_pfq([1,1],[],CC(-1)),
                 CC("[0.59634736232319407434107849937 +/- 2.01e-30]"))

  # 0!^2 - 1!^2 + 2!^2 - 3!^2 + ... = 0.66...
  @test overlaps(hypergeometric_pfq([1,1,1],[],CC(-1)),
                 CC("[0.66809132637777776543301498750 +/- 1.58e-30]"))

  # 0!^3 - 1!^3 + 2!^3 - 3!^3 + ... = 0.72...
  @test overlaps(hypergeometric_pfq([1,1,1,1],[],CC(-1)),
                 CC("[0.72374994418155305222032415740 +/- 5.33e-30]"))
end

end

function run_bench(a, b, z, p, fx = "", fy = "")

  A = map(QQ, a)
  B = map(QQ, b)
  println("---- $A $B $z digits=$p ----")

  p = trunc(Int, p*3.322)
  zz = nacb(narb(QQ(real(z)), p), narb(QQ(imag(z)), p))
  @time f = compute_pfq(A, B, zz, p)
  @assert isempty(fx) || overlaps(set!(narb(), fx, p), real(f))
  @assert isempty(fy) || overlaps(set!(narb(), fy, p), imag(f))
  bits_lost = p - precision(f)
  return bits_lost
end

function run_benchmarks()
println()
println("************************** benchmarks *******************************")

  t1 = time()
  bits_lost = 0

  bits_lost += run_bench([1//2,1//3,1//5,1//4],[8//7,1//6,1//9], 1//10, 20, "[1.04165510107809208771 +/- 2.43e-21]", "[+/- 2.07e-24]")

  bits_lost += run_bench([1//2,1//3,1//5,1//4],[8//7,1//6,1//9],-1+im//300, 200, "[0.7361464107568699653654181932021195729016525520251669994612965514688945567873978145804140256340981534032171785111170341312718204804372620880207341602115123749944939921894542469460886808720225689828743232125 +/- 7.05e-206]", "[0.0006091142434506593619507708486086804762150980760156332916347162351588534898343587334997033663903136796976875663495423707271793534625090625898822976732269967397826684043934917897380440587935220808422516 +/- 1.57e-203]")
  bits_lost += run_bench([1//2,1//3,1//5,1//4],[8//7,1//6,1//9],1//300+im, 200, "[0.853498211060579954490363364784481268197877903709746131856530610099120791429240016786586006359060701719614321941609976993732551828969874878632105881859561180925930687147324142274043333887472085663429960 +/- 4.43e-202]", "[0.302235718305177081880438308463901805943967344220865402810895162185350037204941683860825775142412188198314838346748247659487319633322537866130218549109342613855275799727201394499484865564020106350734127 +/- 2.99e-202]")

  bits_lost += run_bench([1//2,1//3,1//5,1//4],[8//7,1//6,1//9], 11//10, 100, "[2.09473837006811761028534487482955095697725399923738208299473100086979606619943028860894705913795820 +/- 2.23e-99]", "[-1.32063424488884879145689203640820495365830310353133670921975318942279074195339938129952186569083161 +/- 3.11e-99]")
  bits_lost += run_bench([1//2,1//3,1//5,1//4],[8//7,1//6,409//420], 11//10, 100, "[1.10012992690145570229907705218466524693847390580609878681247630305435529779555060922581931979395195 +/- 7.02e-100]", "[-0.019650631936840220248345681395791981619880186753148199866219643526635725512689324474868785073820769 +/- 1.99e-100]")

  bits_lost += run_bench([1//2,1//3,1//5,1//4],[8//7,1//6,1//9], 21//10+im, 100, "[0.71474921490254889156840981949000811062050067188575315042564339656971171954121748902674586037861411 +/- 4.77e-99]", "[0.9781249447537276082866697472831284801654120848168295066054521588403569889034210483850430568298298 +/- 3.56e-98]")

  bits_lost += run_bench([1//2,1//3,1//5,1//4],[8//7,1//6,10//9], 21//10, 100, "[1.09611960714180653165793929034538750471702782740589204327790995048489804128469673421588286818944088 +/- 2.55e-99]", "[-0.1103381390329406104849629727151339938760273060221215999605687131539917991045753083253694103910355 +/- 4.15e-98]")

  bits_lost += run_bench([1//2,1//3,3//2,3//2],[8//7,1//6,409//420], 21//10*im, 100, "[-0.3836936156357715109532470008774342030087022878339141902193492463878637848072764973662114392423830 +/- 3.82e-98]", "[0.3510560688879407205546122233658260321782263643876690384980560383966660407148203677984173909952239 +/- 4.45e-98]")

  bits_lost += run_bench([1//2,1//3,3//2,3//2],[8//7,1//6,7//2+1//4], 21//10*im, 100, "[0.498286277362499154454672287676933977606853722206194319858084629633516142945046940166960257970613 +/- 3.22e-97]", "[0.5388483851104630403622090287525048639997504587807914458817597010957463072020568473516818023683285 +/- 4.40e-98]")

  bits_lost += run_bench([1//2,1//3,3//2,3//2],[8//7,1//6,409//420], 21//10, 100, "[-0.3778888881782393774553127309816076760288807629794011361623259687233181125952890271460921393790923 +/- 1.00e-98]", "[1.0324057197731830866862412722661912510450592489443071095814531303361891341977524868140295288542887 +/- 3.18e-98]")

  bits_lost += run_bench([1//2,1//3,3//2,3//2],[8//7,1//6,7//2+1//4], 21//10, 100, "[0.5587476752887922757498928650974538228321146670224016380361419696854589138256157940556936141300072 +/- 7.32e-98]", "[-2.649787031624682836178763161718415233309248545654569640464561282640046441505973680884328801960931 +/- 1.55e-97]")
  bits_lost += run_bench([1//2,1//3,3//2,5//2],[8//7,1//6,7//2+1//4], 21//10, 100, "[-1.799478585201052861552266948753576293798081804909492183540383822249557260558139389282027739306081 +/- 2.29e-97]", "[-2.0550213095299542518482540717224091776461005489751204129826789474740506081260337268088507033546378 +/- 5.82e-98]")
  bits_lost += run_bench([1//2,1//3,1//3,3//2,3//2],[8//7,1//6,7//2+1//5,1//4], 21//10, 100, "[-0.116857050064488324880005256242593338635544257506370503116752352065017363357264920954229325317036 +/- 5.51e-97]", "[-3.515406819990606948896449062181976165047281121147375761111710596644569438443262697365115597303044 +/- 2.31e-97]")
  bits_lost += run_bench([1//2,1//3,1//3,3//2,3//2,4//3],[8//7,1//6,7//2,1//4,1//5], 21//10, 100, "[-19.9404485154538633201163293612929060982693170333095030435564215171671004250390889987049810825485 +/- 2.54e-95]", "[17.0022670479727611109224091973883810218739134352957349354752346919984236719406436188027263944349 +/- 4.04e-95]")
  bits_lost += run_bench([1//2,1//3,1//3,3//2,3//2,4//3],[8//7,1//6,7//2,1//4,5+1//5], 101//100, 100)
  bits_lost += run_bench([1,1,2,2,3],[3//2,3//2,3//2,3//2], 21//10, 100, "[-0.51585621412817788953636435138740344792545648221588737349075642969025142618998460938638704600278106 +/- 3.91e-99]", "[-0.05470891931824783630344691907005460888039986744824234445189809615982014632842245019368879129590597 +/- 3.12e-99]")

  bits_lost += run_bench([-2,1,2,2,3],[4//3,4//3,4//3,4//3], 21//10, 100, "[25.724630102040816326530612244897959183673469387755102040816326530612244897959183673469387755102041 +/- 2.05e-97]")

  println("total: ", time()-t1)
  println("bits lost: ", bits_lost)
end

function run_heuristics()
println()
println("************************** heuristics *******************************")
  for (z, limit) in [(1, 50000), (1-1//50, 11000)]
    d = 10
    while (d = 2*d) <= limit
      m = trunc(Int, d*4.0)
      N = trunc(Int, d*1.0)
      p = trunc(Int, d*3.322)
      println()
      println("--------- $d digits, z = $z ----------")

      @time s1 = heuristic_hyp(map(QQ,[1//2,1//3,1//4,3//5]),
                               map(QQ, [11//7,7//6,1//8]), QQ(z),
                               m, N, p)

      @time s2 = hypergeometric_pfq([1//2,1//3,1//4,3//5],
                                    [11//7,7//6,1//8], AcbField(p)(z))
      s2 = nacb(s2)
      @show sub(s1, s2, p+60)
      if !overlaps(s1, s2)
        println("oops: ", sub(s2, s1, p+100))
        error("bad answer")
      end
    end
  end
end

run_tests()
run_benchmarks()
run_heuristics()

nothing

