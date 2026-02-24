#' Create a Tuning Grid of Survival Learners
#'
#' Dynamically generates custom wrapper functions for a specified base learner
#' across a grid of hyperparameters.
#'
#' @param base_learner Character string of the base learner function name (e.g., "surv.gbm").
#' @param grid_params List of numeric/character vectors containing hyperparameter values.
#' @return A character vector of the newly generated function names.
#' @export
#' @keywords internal
create_surv_grid <- function(base_learner, grid_params) {

  # 1. Create a data.frame of all possible hyperparameter combinations
  param_grid <- expand.grid(grid_params, stringsAsFactors = FALSE)

  # Initialize a vector to store the names of the new functions we are about to make
  generated_learners <- character(nrow(param_grid))

  for (i in seq_len(nrow(param_grid))) {

    # 2. Create a clean, unique name for this specific learner (e.g., "surv.gbm_1")
    learner_name <- paste0(base_learner, "_", i)
    generated_learners[i] <- learner_name

    # Extract the specific parameters for this exact row as a list
    specific_params <- as.list(param_grid[i, , drop = FALSE])

    # 3. Dynamically construct the new wrapper function
    # We use local() to lock in the specific parameters for this iteration of the loop
    new_fun <- local({

      .base_learner <- base_learner
      .specific_params <- specific_params

      # This is the actual function that will be spawned
      function(time, event, X, newX, new.times, obsWeights = NULL, id = NULL, ...) {

        # Grab the original base algorithm (e.g., surv.gbm)
        base_fun <- get(.base_learner, mode = "function")

        # Set up the standard required arguments
        args_to_pass <- list(
          time = time,
          event = event,
          X = X,
          newX = newX,
          new.times = new.times,
          obsWeights = obsWeights,
          id = id
        )

        # Merge the standard args, the grid's tuned hyperparameters, and any ...
        call_args <- c(args_to_pass, .specific_params, list(...))

        # Execute the base learner with these perfectly locked-in settings
        do.call(base_fun, call_args)
      }
    })

    # 4. Inject the newly created function into the global R environment
    assign(learner_name, new_fun, envir = parent.frame())
  }

  # Return the vector of names so it can be combined into the SuperSurv library list
  return(generated_learners)
}
