module FITSIO

export FITSFile,
       fits_clobber_file,
       fits_close_file,
       fits_create_ascii_tbl,
       fits_create_binary_tbl,
       fits_create_file,
       fits_create_img,
       fits_delete_file,
       fits_delete_key,
       fits_delete_record,
       fits_delete_rows,
       fits_file_mode,
       fits_file_name,
       fits_get_col_repeat,
       fits_get_hdrspace,
       fits_get_hdu_num,
       fits_get_img_size,
       fits_get_num_cols,
       fits_get_num_hdus,
       fits_get_num_rows,
       fits_get_num_rowsll,
       fits_get_rowsize,
       fits_hdr2str,
       fits_insert_rows,
       fits_movabs_hdu,
       fits_movrel_hdu,
       fits_movnam_hdu,
       fits_open_data,
       fits_open_file,
       fits_open_image,
       fits_open_table,
       fits_read_col,
       fits_read_keyn,
       fits_read_keyword,
       fits_read_pix,
       fits_read_record,
       fits_write_col,
       fits_write_key,
       fits_write_pix,
       fits_write_record

export fitsread

# ------------------------------------
# High-level interface
export FITS,
       HDU,
       ImageHDU

import Base: getindex, length, show, read, write, close, ndims, size
# ------------------------------------

using BinDeps
@BinDeps.load_dependencies

type FITSFile
    ptr::Ptr{Void}

    function FITSFile(ptr::Ptr{Void})
        f = new(ptr)
        finalizer(f, close)
        f
    end
end

function fits_assert_open(f::FITSFile)
    if f.ptr == C_NULL
        error("attempt to access closed FITS file")
    end
end

function fits_get_errstatus(status::Int32)
    msg = Array(Uint8, 31)
    ccall((:ffgerr,libcfitsio), Void, (Int32,Ptr{Uint8}), status, msg)
    bytestring(convert(Ptr{Uint8},msg))
end

function fits_assert_ok(status::Int32)
    if status != 0
        error(fits_get_errstatus(status))
    end
end

# constants

_cfitsio_bitpix{T<:Integer}(::Type{T}) = int32(8*sizeof(T))
_cfitsio_bitpix{T<:FloatingPoint}(::Type{T}) = int32(-8*sizeof(T))

_cfitsio_datatype(::Type{Uint8})      = int32(11)
_cfitsio_datatype(::Type{Int8})       = int32(12)
_cfitsio_datatype(::Type{Bool})       = int32(14)
_cfitsio_datatype{T<:String}(::Type{T}) = int32(16)
_cfitsio_datatype(::Type{Uint16})     = int32(20)
_cfitsio_datatype(::Type{Int16})      = int32(21)
_cfitsio_datatype(::Type{Uint32})     = int32(30)
_cfitsio_datatype(::Type{Int32})      = int32(31)
_cfitsio_datatype(::Type{Float32})    = int32(42)
_cfitsio_datatype(::Type{Int64})      = int32(81)
_cfitsio_datatype(::Type{Float64})    = int32(82)
_cfitsio_datatype(::Type{Complex64})  = int32(83)
_cfitsio_datatype(::Type{Complex128}) = int32(163)

function hdu_int_to_type(hdu_type_int)
    if hdu_type_int == 0
        return :image_hdu
    elseif hdu_type_int == 1
        return :ascii_table
    elseif hdu_type_int == 2
        return :binary_table
    end

    :unknown
end

const mode_strs = [int32(0)=>"READONLY", int32(1)=>"READWRITE"]

const bitpix_to_str = [int32(8)=>"BYTE_IMG", int32(16)=>"SHORT_IMG",
                       int32(32)=>"LONG_IMG", int32(64)=>"LONGLONG_IMG",
                       int32(-32)=>"FLOAT_IMG", int32(-64)=>"DOUBLE_IMG"]
const bitpix_to_type = [int32(8)=>Uint8, int32(16)=>Int16,
                        int32(32)=>Int32, int32(64)=>Int64,
                        int32(-32)=>Float32, int32(-64)=>Float64]

# General-purpose functions

for (a,b,T) in ((:fits_file_mode,     "ffflmd",  :Cint),
                (:fits_get_num_cols,  "ffgncl",  :Cint),
                (:fits_get_num_hdus,  "ffthdu",  :Cint),
                (:fits_get_num_rows,  "ffgnrw",  :Clong),
                (:fits_get_num_rowsll,"ffgnrwll",:Clonglong),
                (:fits_get_rowsize,   "ffgrsz",  :Cint))
    @eval begin
        function ($a)(f::FITSFile)
            result = $T[0]
            status = Int32[0]
            ccall(($b,libcfitsio), Int32,
                  (Ptr{Void}, Ptr{$T}, Ptr{Int32}),
                  f.ptr, result, status)
            fits_assert_ok(status[1])
            result[1]
        end
    end
end

# file access

function fits_create_file(filename::String)
    ptr = Array(Ptr{Void}, 1)
    status = Int32[0]
    ccall((:ffinit,libcfitsio),
        Int32, (Ptr{Ptr{Void}},Ptr{Uint8},Ptr{Int32}),
        ptr, bytestring(filename), status)
    fits_assert_ok(status[1])
    FITSFile(ptr[1])
end

fits_clobber_file(filename::String) = fits_create_file("!"*filename)

for (a,b) in ((:fits_open_data, "ffdopn"),
              (:fits_open_file, "ffopen"),
              (:fits_open_image,"ffiopn"),
              (:fits_open_table,"fftopn"))
    @eval begin
        function ($a)(filename::String)
            ptr = Array(Ptr{Void}, 1)
            mode = int32(0) # readonly
            status = Int32[0]
            ccall(($b,libcfitsio), Int32,
                  (Ptr{Ptr{Void}},Ptr{Uint8},Int32,Ptr{Int32}),
                  ptr, bytestring(filename), mode, status)
            fits_assert_ok(status[1])
            FITSFile(ptr[1])
        end
    end
end

for (a,b) in ((:fits_close_file, "ffclos"),
              (:fits_delete_file,"ffdelt"))
    @eval begin
        function ($a)(f::FITSFile)
            
            # fits_close_file is called during garbage collection, but the
            # file may already be closed.
            if f.ptr == C_NULL
                return
            end

            status = Int32[0]
            ccall(($b,libcfitsio), Int32,
                  (Ptr{Void},Ptr{Int32}),
                  f.ptr, status)
            fits_assert_ok(status[1])
            f.ptr = C_NULL
            nothing
        end
    end
end

close(f::FITSFile) = fits_close_file(f)

function fits_file_name(f::FITSFile)
    value = Array(Uint8, 1025)
    status = Int32[0]
    ccall((:ffflnm,libcfitsio), Int32,
          (Ptr{Void},Ptr{Uint8},Ptr{Int32}),
          f.ptr, value, status)
    fits_assert_ok(status[1])
    bytestring(convert(Ptr{Uint8}, value))
end

# header keywords

function fits_get_hdrspace(f::FITSFile)
    keysexist = Int32[0]
    morekeys = Int32[0]
    status = Int32[0]
    ccall((:ffghsp,libcfitsio), Int32,
        (Ptr{Void},Ptr{Int32},Ptr{Int32},Ptr{Int32}),
        f.ptr, keysexist, morekeys, status)
    fits_assert_ok(status[1])
    (keysexist[1], morekeys[1])
end

function fits_read_keyword(f::FITSFile, keyname::String)
    value = Array(Uint8, 71)
    comment = Array(Uint8, 71)
    status = Int32[0]
    ccall((:ffgkey,libcfitsio), Int32,
        (Ptr{Void},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Int32}),
        f.ptr, bytestring(keyname), value, comment, status)
    fits_assert_ok(status[1])
    bytestring(convert(Ptr{Uint8},value)), bytestring(convert(Ptr{Uint8},comment))
end

function fits_read_record(f::FITSFile, keynum::Int)
    card = Array(Uint8, 81)
    status = Int32[0]
    ccall((:ffgrec,libcfitsio), Int32,
        (Ptr{Void},Int32,Ptr{Uint8},Ptr{Int32}),
        f.ptr, keynum, card, status)
    fits_assert_ok(status[1])
    bytestring(convert(Ptr{Uint8},card))
end

function fits_read_keyn(f::FITSFile, keynum::Int)
    keyname = Array(Uint8, 9)
    value = Array(Uint8, 71)
    comment = Array(Uint8, 71)
    status = Int32[0]
    ccall((:ffgkyn,libcfitsio), Int32,
        (Ptr{Void},Int32,Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Int32}),
        f.ptr, keynum, keyname, value, comment, status)
    fits_assert_ok(status[1])
    bytestring(convert(Ptr{Uint8},keyname)),
    bytestring(convert(Ptr{Uint8},value)),
    bytestring(convert(Ptr{Uint8},comment))
end

function fits_write_key(f::FITSFile, keyname::String, value::Union(Number,String), comment::String)
    cvalue = isa(value,String) ?  bytestring(value) :
             isa(value,Bool) ? [int32(value)] : [value]
    status = Int32[0]
    ccall((:ffpky,libcfitsio), Int32,
        (Ptr{Void},Int32,Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Int32}),
        f.ptr, _cfitsio_datatype(typeof(value)), bytestring(keyname),
        cvalue, bytestring(comment), status)
    fits_assert_ok(status[1])
end

function fits_write_record(f::FITSFile, card::String)
    status = Int32[0]
    ccall((:ffprec,libcfitsio), Int32,
        (Ptr{Void},Ptr{Uint8},Ptr{Int32}),
        f.ptr, bytestring(card), status)
    fits_assert_ok(status[1])
end

function fits_delete_record(f::FITSFile, keynum::Int)
    status = Int32[0]
    ccall((:ffdrec,libcfitsio), Int32,
        (Ptr{Void},Int32,Ptr{Int32}),
        f.ptr, keynum, status)
    fits_assert_ok(status[1])
end

function fits_delete_key(f::FITSFile, keyname::String)
    status = Int32[0]
    ccall((:ffdkey,libcfitsio), Int32,
        (Ptr{Void},Ptr{Uint8},Ptr{Int32}),
        f.ptr, bytestring(keyname), status)
    fits_assert_ok(status[1])
end

function fits_hdr2str(f::FITSFile, nocomments::Bool=false)
    status = Int32[0]
    header = Array(Ptr{Uint8}, 1)
    nkeys = Int32[0]
    ccall((:ffhdr2str, libcfitsio), Ptr{Ptr{Uint8}},
          (Ptr{Void},Int32, Ptr{Ptr{Uint8}}, Int32,
           Ptr{Ptr{Uint8}}, Ptr{Int32}, Ptr{Int32}),
          f.ptr, nocomments, &C_NULL, 0, header, nkeys, status)
    result = bytestring(header[1])

    # free header pointer allocated by cfitsio (result is a copy)
    ccall((:fffree, libcfitsio), Ptr{Int32},
          (Ptr{Uint8}, Ptr{Int32}),
          header[1], status)
    fits_assert_ok(status[1])
    result
end

# HDU functions

for (a,b) in ((:fits_movabs_hdu,"ffmahd"),
              (:fits_movrel_hdu,"ffmrhd"))
    @eval begin
        function ($a)(f::FITSFile, hduNum::Integer)
            hdu_type = Int32[0]
            status = Int32[0]
            ccall(($b,libcfitsio), Int32,
                  (Ptr{Void}, Int32, Ptr{Int32}, Ptr{Int32}),
                  f.ptr, hduNum, hdu_type, status)
            fits_assert_ok(status[1])
            hdu_int_to_type(hdu_type[1])
        end
    end
end

function fits_movnam_hdu(f::FITSFile, extname::String, extver::Integer=0,
                         hdu_type::Integer=-1)
    status = Int32[0]
    ccall((:ffmnhd,libcfitsio), Int32,
          (Ptr{Void}, Int32, Ptr{Uint8}, Int32, Ptr{Int32}),
          f.ptr, int32(hdu_type), bytestring(extname), int32(extver), status)
    fits_assert_ok(status[1])
end

function fits_get_hdu_num(f::FITSFile)
    hdunum = Int32[0]
    ccall((:ffghdn,libcfitsio), Int32,
          (Ptr{Void},Ptr{Int32}),
          f.ptr, hdunum)
    hdunum[1]
end

function fits_get_hdu_type(f::FITSFile)
    hdutype = Int32[0]
    status = Int32[0]
    ccall((:ffghdt, libcfitsio), Int32,
          (Ptr{Void}, Ptr{Int32}, Ptr{Int32}),
          f.ptr, hdutype, status)
    fits_assert_ok(status[1])
    hdu_int_to_type(hdutype[1])
end

# primary array or IMAGE extension

function fits_get_img_type(f::FITSFile)
    bitpix = Int32[0]
    status = Int32[0]
    ccall((:ffgidt,libcfitsio), Int32,
        (Ptr{Void},Ptr{Int32},Ptr{Int32}),
        f.ptr, bitpix, status)
    fits_assert_ok(status[1])
    bitpix[1]
end

function fits_get_img_dim(f::FITSFile)
    ndim = Int32[0]
    status = Int32[0]
    ccall((:ffgidm,libcfitsio), Int32,
        (Ptr{Void},Ptr{Int32},Ptr{Int32}),
        f.ptr, ndim, status)
    fits_assert_ok(status[1])
    ndim[1]
end

function fits_get_img_size(f::FITSFile)
    ndim = fits_get_img_dim(f)
    naxes = zeros(Int, ndim)
    status = Int32[0]
    ccall((:ffgisz,libcfitsio), Int32,
        (Ptr{Void},Int32,Ptr{Int},Ptr{Int32}),
        f.ptr, ndim, naxes, status)
    fits_assert_ok(status[1])
    naxes
end

function fits_create_img(f::FITSFile, t::Type, naxes::Vector{Int})
    status = Int32[0]
    ccall((:ffcrim,libcfitsio), Int32,
        (Ptr{Void},Int32,Int32,Ptr{Int},Ptr{Int32}),
        f.ptr, _cfitsio_bitpix(t), length(naxes), naxes, status)
    fits_assert_ok(status[1])
end

function fits_write_pix{T}(f::FITSFile, fpixel::Vector{Int}, nelements::Int, data::Array{T})
    status = Int32[0]
    ccall((:ffppx,libcfitsio), Int32,
        (Ptr{Void},Int32,Ptr{Int},Int,Ptr{Void},Ptr{Int32}),
        f.ptr, _cfitsio_datatype(T), fpixel, nelements, data, status)
    fits_assert_ok(status[1])
end
fits_write_pix(f::FITSFile, data::Array) = fits_write_pix(f, ones(Int,length(size(data))), length(data), data)

function fits_read_pix{T}(f::FITSFile, fpixel::Vector{Int}, nelements::Int, nullval::T, data::Array{T})
    anynull = Int32[0]
    status = Int32[0]
    ccall((:ffgpxv,libcfitsio), Int32,
        (Ptr{Void},Int32,Ptr{Int},Int,Ptr{Void},Ptr{Void},Ptr{Int32},Ptr{Int32}),
        f.ptr, _cfitsio_datatype(T), fpixel, nelements, &nullval, data, anynull, status)
    fits_assert_ok(status[1])
    anynull[1]
end

function fits_read_pix{T}(f::FITSFile, fpixel::Vector{Int}, nelements::Int, data::Array{T})
    anynull = Int32[0]
    status = Int32[0]
    ccall((:ffgpxv,libcfitsio), Int32,
        (Ptr{Void},Int32,Ptr{Int},Int,Ptr{Void},Ptr{Void},Ptr{Int32},Ptr{Int32}),
        f.ptr, _cfitsio_datatype(T), fpixel, nelements, C_NULL, data, anynull, status)
    fits_assert_ok(status[1])
    anynull[1]
end
fits_read_pix(f::FITSFile, data::Array) = fits_read_pix(f, ones(Int,length(size(data))), length(data), data)

function fits_read_subset{T}(f::FITSFile,
                             fpixel::Vector{Clong},
                             lpixel::Vector{Clong},
                             inc::Vector{Clong},
                             data::Array{T})
    anynull = Cint[0]
    status = Cint[0]
    ccall((:ffgsv,libcfitsio), Cint,
          (Ptr{Void},Cint,Ptr{Clong},Ptr{Clong},Ptr{Clong},Ptr{Void},Ptr{Void},
           Ptr{Cint},Ptr{Cint}),
          f.ptr, _cfitsio_datatype(T), fpixel, lpixel, inc, C_NULL, data,
          anynull, status)
    fits_assert_ok(status[1])
    anynull[1]
end

function fitsread(filename::String)
    f = fits_open_file(filename)
    s = fits_get_img_size(f)
    a = Array(Float64, s...)
    fits_read_pix(f, a)
    fits_close_file(f)
    a'
end

# ASCII/binary tables

# The three fields are: ttype, tform, tunit (CFITSIO's terminology)
typealias ColumnDef (ASCIIString, ASCIIString, ASCIIString)

for (a,b) in ((:fits_create_binary_tbl, 2),
              (:fits_create_ascii_tbl,  1))
    @eval begin
        function ($a)(f::FITSFile, numrows::Integer,
                      coldefs::Array{ColumnDef}, extname::String)

            ntype = length(coldefs)
            ttype = map((x) -> pointer(x[1].data), coldefs)
            tform = map((x) -> pointer(x[2].data), coldefs)
            tunit = map((x) -> pointer(x[3].data), coldefs)
            status = Int32[0]

            ccall(("ffcrtb", libcfitsio), Int32,
                  (Ptr{Void}, Int32, Int64, Int32,
                   Ptr{Ptr{Uint8}}, Ptr{Ptr{Uint8}},
                   Ptr{Ptr{Uint8}}, Ptr{Uint8}, Ptr{Int32}),
                  f.ptr, $b, convert(Int32, numrows), ntype,
                  ttype, tform, tunit, bytestring(extname),
                  status)
            fits_assert_ok(status[1])
        end
    end
end

function fits_get_col_repeat(f::FITSFile, colnum::Integer)
    typecode = Int32[0]
    repeat = Int64[0]
    width = Int64[0]
    status = Int32[0]

    ccall((:ffgtclll, libcfitsio), Int32,
          (Ptr{Void}, Int32, Ptr{Int32}, Ptr{Int64}, Ptr{Int64}, Ptr{Int32}),
          f.ptr, convert(Int32, colnum), typecode, repeat, width, status)

    (repeat[1], width[1])
end

function fits_read_col{T}(f::FITSFile,
                          ::Type{T},
                          colnum::Integer,
                          firstrow::Integer,
                          firstelem::Integer,
                          data::Array{T})

    anynull = Int32[0]
    status = Int32[0]

    firstrow64 = convert(Int64, firstrow)
    firstelem64 = convert(Int64, firstelem)
    nelements64 = convert(Int64, length(data))

    if isa(T, Type{String}) || isa(T, Type{ASCIIString})

        # Make sure there is enough room for each string
        repcount, width = fits_get_col_repeat(f, colnum)
        for i in 1:length(data)
            # We need to call `repeat' N times in order to allocate N
            # strings
            data[i] = repeat(" ", repcount)
        end

        ccall((:ffgcvs, libcfitsio), Int32,
              (Ptr{Void}, Int32, Int64, Int64, Int64,
               Ptr{Uint8}, Ptr{Ptr{Uint8}}, Ptr{Int32}, Ptr{Int32}),
              f.ptr, convert(Int32, colnum),
              firstrow64, firstelem64, nelements64,
              "", data, anynull, status)

        # Truncate the strings to the first NULL character (if present)
        for idx in 1:length(data)
            zeropos = search(data[idx], '\0')
            if zeropos >= 1
                data[idx] = (data[idx])[1:(zeropos-1)]
            end
        end

    else

        ccall((:ffgcv,libcfitsio), Int32,
              (Ptr{Void}, Int32, Int32, Int64, Int64, Int64,
               Ptr{T}, Ptr{T}, Ptr{Int32}, Ptr{Int32}),
              f.ptr, _cfitsio_datatype(T), convert(Int32, colnum),
              firstrow64, firstelem64, nelements64,
              T[0], data, anynull, status)

    end

    fits_assert_ok(status[1])

end

function fits_write_col{T}(f::FITSFile,
                           ::Type{T},
                           colnum::Integer,
                           firstrow::Integer,
                           firstelem::Integer,
                           data::Array{T})

    firstrow64 = convert(Int64, firstrow)
    firstelem64 = convert(Int64, firstelem)
    nelements64 = convert(Int64, length(data))
    status = Int32[0]

    if  isa(T, Type{String}) || isa(T, Type{ASCIIString})

        ccall((:ffpcls, libcfitsio), Int32,
              (Ptr{Void}, Int32, Int64, Int64, Int64,
               Ptr{Ptr{Uint8}}, Ptr{Int32}),
              f.ptr, convert(Int32, colnum),
              firstrow64, firstelem64, nelements64,
              data, status)

    else

        ccall((:ffpcl,libcfitsio), Int32,
              (Ptr{Void}, Int32, Int32, Int64, Int64, Int64,
               Ptr{T}, Ptr{Int32}),
              f.ptr, _cfitsio_datatype(T), convert(Int32, colnum),
              firstrow64, firstelem64, nelements64,
              data, status)

    end

    fits_assert_ok(status[1])

end

for (a,b) in ((:fits_insert_rows, "ffirow"),
              (:fits_delete_rows, "ffdrow"))
    @eval begin
        function ($a)(f::FITSFile, firstrow::Integer, nrows::Integer)
            firstrow64 = convert(Int64, firstrow)
            nrows64 = convert(Int64, nrows)
            status = Int32[0]

            ccall(($b,libcfitsio), Int32,
                  (Ptr{Void}, Int64, Int64, Ptr{Int32}),
                  f.ptr, firstrow64, nrows64, status)
            fits_assert_ok(status[1])
        end
    end
end

# -----------------------------------------------------------------------------
# High-level interface: Types

abstract HDU

# FITS is analagous to FITSFile, but holds a reference to all of its
# HDU objects. This is so that only a single HDU object is created for
# each extension in the file. It also allows a FITS object to tell
# previously created HDUs about events that happen to the file, such 
# as deleting extensions. This could be done by, e.g., setting ext=-1 in
# the HDU object.

type FITS
    fitsfile::FITSFile
    filename::String
    mode::String
    hdus::Dict{Int, HDU}

    function FITS(filename::String, mode::String="r")
        f = (mode == "r" ? fits_open_file(filename) :
             mode == "r+" && isfile(filename) ? error("r/w mode not yet implemented"):
             mode == "r+" ? fits_create_file(filename) :
             mode == "w" ? fits_create_file("!"*filename) :
             error("invalid open mode: $mode"))

        new(f, filename, mode, Dict{Int, HDU}())
    end
end

# TODO : Cache metadata such as extname, extver, image size and data type?
#        This might allow faster access for size(ImageHDU) and ndim(ImageHDU)
type ImageHDU <: HDU
    fitsfile::FITSFile
    ext::Int
end

type TableHDU <: HDU
    fitsfile::FITSFile
    ext::Int
end

type AsciiHDU <: HDU
    fitsfile::FITSFile
    ext::Int
end

# -----------------------------------------------------------------------------
# FITS

function length(f::FITS)
    fits_assert_open(f.fitsfile)
    fits_get_num_hdus(f.fitsfile)
end

endof(f::FITS) = length(f)

function show(io::IO, f::FITS)
    fits_assert_open(f.fitsfile)

    print(io, "file: ", f.filename, "\n")
    print(io, "mode: ", f.mode, "\n")
    print(io, "extnum exttype         extname\n")

    for i = 1:length(f)
        hdutype = fits_movabs_hdu(f.fitsfile, i)
        extname = ""
        try
            extname = fits_read_keyword(f.fitsfile, "EXTNAME")[1]
        catch
            try
                extname = fits_read_keyword(f.fitsfile, "HDUNAME")[1]
            catch
            end
        end
        @printf io "%-6d %-15s %s\n" i hdutype extname
    end
end

# Returns HDU object based on extension number
function getindex(f::FITS, i::Int)
    fits_assert_open(f.fitsfile)

    if haskey(f.hdus, i)
        return f.hdus[i]
    end

    hdutype = fits_movabs_hdu(f.fitsfile, i)
    f.hdus[i] = (hdutype == :image_hdu ? ImageHDU(f.fitsfile, i) :
                 hdutype == :binary_table ? TableHDU(f.fitsfile, i) :
                 hdutype == :ascii_table ? AsciiHDU(f.fitsfile, i) :
                 error("bad HDU type"))
    return f.hdus[i]
end

# Returns HDU based on hduname, version
function getindex(f::FITS, name::String, ver::Int=0)
    fits_assert_open(f.fitsfile)
    fits_movnam_hdu(f.fitsfile, name, ver)
    i = fits_get_hdu_num(f.fitsfile)

    if haskey(f.hdus, i)
        return f.hdus[i]
    end

    hdutype = fits_get_hdu_type(f.fitsfile)
    f.hdus[i] = (hdutype == :image_hdu ? ImageHDU(f.fitsfile, i) :
                 hdutype == :binary_table ? TableHDU(f.fitsfile, i) :
                 hdutype == :ascii_table ? AsciiHDU(f.fitsfile, i) :
                 error("bad HDU type"))
    return f.hdus[i]
end


function close(f::FITS)
    fits_assert_open(f.fitsfile)
    fits_close_file(f.fitsfile)
    f.filename = ""
    f.mode = ""
    empty!(f.hdus)
    nothing
end

# -----------------------------------------------------------------------------
# ImageHDU methods

# Display the image datatype and dimensions
function show(io::IO, hdu::ImageHDU)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    bitpix = fits_get_img_type(hdu.fitsfile)
    sz = fits_get_img_size(hdu.fitsfile)
    @printf io "file: %s\nextension: %d\ntype: IMAGE\nimage info:\n  bitpix: %d\n  size: %s" fits_file_name(hdu.fitsfile) hdu.ext bitpix tuple(sz...)
end

# Get image dimensions
function ndims(hdu::ImageHDU)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    fits_get_img_dim(hdu.fitsfile)
end

# Get image size
function size(hdu::ImageHDU)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    sz = fits_get_img_size(hdu.fitsfile)
    tuple(sz...)
end

function size(hdu::ImageHDU, i::Integer)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    sz = fits_get_img_size(hdu.fitsfile)
    sz[i]
end

# Read a full image from an HDU
# TODO: Correct support for BSCALE'd images
function read(hdu::ImageHDU)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    sz = fits_get_img_size(hdu.fitsfile)
    bitpix = fits_get_img_type(hdu.fitsfile)
    data = Array(bitpix_to_type[bitpix], sz...)
    fits_read_pix(hdu.fitsfile, data)
    data
end

# Read all or a subset of an HDU
function getindex(hdu::ImageHDU, rs::Range...)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)

    # construct first, last and step vectors
    firsts = Clong[first(r) for r in rs]
    lasts = Clong[last(r) for r in rs]
    steps = Clong[step(r) for r in rs]

    # construct output array
    bitpix = fits_get_img_type(hdu.fitsfile)
    
    datasz = [length(r) for r in rs]
    data = Array(bitpix_to_type[bitpix], datasz...)
    fits_read_subset(hdu.fitsfile, firsts, lasts, steps, data)
    data
end

# Add a new ImageHDU to a FITS object
function write{T}(f::FITS, data::Array{T})
    fits_assert_open(f.fitsfile)
    s = size(data)
    fits_create_img(f.fitsfile, T, [s...])
    fits_write_pix(f.fitsfile, ones(Int, length(s)), length(data), data)
    nothing
end

# -----------------------------------------------------------------------------

end # module
