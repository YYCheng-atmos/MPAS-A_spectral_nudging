program pio2
    use pio
    integer, parameter :: MPAS_IO_OFFSET_KIND = PIO_OFFSET_KIND
    integer, parameter :: MPAS_INT_FILLVAL = PIO_FILL_INT
    type (file_desc_t) :: pio_file
    type (Var_desc_t) :: field_desc
    integer (kind=MPAS_IO_OFFSET_KIND) :: frame_number
    call PIO_setframe(pio_file, field_desc, frame_number)
end program
