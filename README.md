# eirene--
A work-in-progress C++ port of the core persistence algorithm from Eirene.jl                                                                                                                                        
                                                                                                                                                                                                                    
Why not use one of the existing C++ persistence packages? Because all of them either don't take general cellular complexes as input or do not have methods for extracting explicit cohomology class representatives.
                                                                                                                                                                                                                    
Does NOT implement the Vietoris Rips -> Cellular Complex algorithm from Eirene.jl, since as I said above I'm only interested in general cell complexes.                                                             
                                                                                                                                                                                                                    
This is still rather untested and should be assumed to be incorrect, and the interface still needs work and documentation, but if you're also interested in the Eirene algorithm it might be helpful.
