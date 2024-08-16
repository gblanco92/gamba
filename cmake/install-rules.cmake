install(
    TARGETS GamBa_exe
    RUNTIME COMPONENT GamBa_Runtime
)

if(PROJECT_IS_TOP_LEVEL)
  include(CPack)
endif()
