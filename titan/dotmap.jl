import Base.getproperty
import Base.setproperty!

mutable struct DotMap
    __dict__
    DotMap() = new(Dict())
end

function DotMap(d::Dict)
   dm = DotMap()
   for (k, v) in d
      dm.__dict__[Symbol(k)] = DotMap(v)
   end
   return dm
end

DotMap(d::Any) = d

# make dots work
function Base.getproperty(obj::DotMap, name::Symbol)
   if name == :__dict__
      return getfield(obj, name)
   else
      return obj.__dict__[name]
   end
end

function Base.setproperty!(obj::DotMap, name::Symbol, x)
   obj.__dict__[name] = x
end

# make dictionary indexing work
function Base.getindex(obj::DotMap, name::Symbol)
   return obj.__dict__[name]
end
Base.getindex(obj::DotMap, name::Any) = Base.getindex(obj, Symbol(name))

function Base.setindex!(obj::DotMap, name::Symbol, x)
   obj.__dict__[name] = x
end

Base.setindex!(obj::DotMap, name, x) = Base.setindex!(obj, Symbol(name), x)

# make iteration work
Base.iterate(obj::DotMap) = Base.iterate(obj.__dict__)
Base.keys(obj::DotMap) = Base.keys(obj.__dict__)
Base.values(obj::DotMap) = Base.values(obj.__dict__)
Base.collect(obj::DotMap) = Base.collect(obj.__dict__)
