".First.lib" <-
function (lib, pkg) 
{
    library.dynam("BsMD", pkg, lib)
    invisible(NULL)
}
