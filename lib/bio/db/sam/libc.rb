require 'rubygems'
require'ffi'
#Binding functions of libc required to do the low level memory managment. 
module LibC
  extend FFI::Library
  ffi_lib FFI::Library::LIBC
  attach_function :free, [ :pointer ], :void
  # call #attach_function to attach to malloc, free, memcpy, bcopy, etc.
end