# Usage: julia bil.jl [argv]
# ======
# This program executes similarly to bin/bil-X.Y.Z.exe
# The implementation is similar to main.c in the C version of bil

macro LIBBIL(); return :("libbil-2.8.8-Debug"); end;

macro Entry_Create(A,B);
  return :(ccall((:Entry_Create, @LIBBIL), Ptr{Cvoid}, (Int,Ptr{Ptr{UInt8}}), $A, $B));
end;

macro Entry_Execute(A);
  return :(ccall((:Entry_Execute, @LIBBIL), Cvoid, (Ref{Cvoid},), $A));
end;

macro Entry_Delete(A);
  return :(ccall((:Entry_Delete, @LIBBIL), Cvoid, (Ref{Cvoid},), $A));
end;

argv = [PROGRAM_FILE]
append!(argv,ARGS)
argc = length(argv)

#entry = @Entry_Create(argc, argv)

@Entry_Execute(@Entry_Create(argc, argv))

#@Entry_Delete(entry)
