
add_executable(${PROJECT_NAME}
  synthese.cpp
# $<TARGET_OBJECTS:${PROJECT_NAME}>
  )

TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${SDL2_LIBRARIES})

#find_library(SDL2 sdl2 sdl2)
#target_link_libraries(${PROJECT_NAME}
#  ${SDL2_LIBRARIES}
#  )

