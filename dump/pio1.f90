program pio1
    use pio
    use pionfatt_mod
    integer, parameter :: MPAS_IO_OFFSET_KIND = PIO_OFFSET
    integer, parameter :: MPAS_INT_FILLVAL = NF_FILL_INT
    type (Var_desc_t) :: field_desc
    integer (kind=MPAS_IO_OFFSET_KIND) :: frame_number
    call PIO_setframe(field_desc, frame_number)
end program
