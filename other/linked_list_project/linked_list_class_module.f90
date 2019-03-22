module linked_list_class_module
    use iso_c_binding
    implicit none

    type :: linked_list_class
        type(c_ptr) :: data
        class(linked_list_class), pointer :: next
    end type

contains
    subroutine linked_list_insert(ptr, data, get_idx)
        type(linked_list_class), pointer, intent(inout) :: ptr
        type(c_ptr), intent(in) :: data
        interface
            function get_idx(data) result(r)
                import
                type(c_ptr), intent(in) :: data
                integer(c_int) :: r
            end function get_idx
        end interface
        type(linked_list_class), pointer, save :: head, tail
        type(linked_list_class), pointer :: tmp, ptr1, ptr2
        integer(c_int) :: idx, head_idx, tail_idx, val1, val2
        logical(c_bool) :: inserted

        if (.not.associated(ptr)) then
            allocate(tmp)
            tmp%data = data
            nullify(tmp%next)
            ptr => tmp
            head => tmp
            tail => tmp
        else
            idx = get_idx(data)
            head_idx = get_idx(head%data)
            tail_idx = get_idx(tail%data)
            if ((idx.ne.head_idx).and.(idx.ne.tail_idx)) then
                if (idx.lt.head_idx) then
                    allocate(tmp)
                    tmp%data = data
                    tmp%next => ptr
                    ptr => tmp
                    head => tmp
                else if (idx.gt.tail_idx) then
                    allocate(tmp)
                    tmp%data = data
                    nullify(tmp%next)
                    tail%next => tmp
                    tail => tmp
                else
                    ptr1 => ptr
                    inserted = .false.
                    do while ((associated(ptr1)).and.(.not.inserted))
                        ptr2 => ptr1%next
                        val1 = get_idx(ptr1%data)
                        val2 = get_idx(ptr2%data)
                        if ((idx.ne.val1).and.(idx.ne.val2)) then
                            if ((idx.gt.val1).and.(idx.lt.val2)) then
                                allocate(tmp)
                                tmp%data = data
                                tmp%next => ptr2
                                ptr1%next => tmp
                                inserted = .true.
                            end if
                        else
                            inserted = .true.
                        end if
                        ptr1 => ptr2
                    end do
                end if
            end if
        end if
    end subroutine linked_list_insert

    subroutine linked_list_free(ptr, free_data)
        type(linked_list_class), intent(inout), pointer :: ptr
        type(c_ptr) :: data
        interface
            subroutine free_data(data)
                import
                type(c_ptr), intent(in) :: data
            end subroutine free_data
        end interface
        class(linked_list_class), pointer :: tmp

        do while (associated(ptr))
            data = ptr%data
            call free_data(data)
            tmp => ptr
            ptr => ptr%next
            deallocate(tmp)
        end do
    end subroutine linked_list_free
end module linked_list_class_module
