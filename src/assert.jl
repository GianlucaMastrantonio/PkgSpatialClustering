
const _USE_ASSERTS = Ref(true)

enable_asserts() = (_USE_ASSERTS[] = true)
disable_asserts() = (_USE_ASSERTS[] = false)
asserts_enabled() = _USE_ASSERTS[]

# ✅ FIXED: use esc to preserve variable scope
macro toggle_assert(cond)
  return quote
    if PkgSpatClust._USE_ASSERTS[]
      @assert $(esc(cond))
    end
  end
end