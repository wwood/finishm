require 'set'
require 'ds'

# Like DS::PriorityQueue except give the ability to define how priority is given
class DS::AnyPriorityQueue < DS::PriorityQueue
  #Create new priority queue. Internaly uses heap to store elements.
  def initialize
    @store = DS::BinaryHeap.new {|parent,child| yield parent.key, child.key}
  end

  def each
    @store.to_a.each do |pair|
      yield pair.value
    end
  end
end



class Bio::AssemblyGraphAlgorithms::BubblyAssembler < Bio::AssemblyGraphAlgorithms::SingleEndedAssembler
  include Bio::FinishM::Logging

  DEFAULT_MAX_BUBBLE_LENGTH = 500
  DEFAULT_BUBBLE_NODE_COUNT_LIMIT = 20 #so, so very 'un-educated' guess
  DEFAULT_BUBBLE_FORK_LIMIT = 20
  DEFAULT_MAX_CYCLES = 1

  def initialize(graph, assembly_options={})
    opts = assembly_options
    opts[:max_bubble_length] ||= DEFAULT_MAX_BUBBLE_LENGTH
    opts[:bubble_node_count_limit] ||= DEFAULT_BUBBLE_NODE_COUNT_LIMIT
    opts[:bubble_fork_limit] ||= DEFAULT_BUBBLE_FORK_LIMIT
    opts[:max_cycles] ||= DEFAULT_MAX_CYCLES
    super graph, opts
  end

  # Starting at a node within a graph, walk through the graph
  # accepting forks, so long as the fork paths converge within some finite
  # length in the graph (the leash length, measured in number of base pairs).
  #
  # Return an Array of Path arrays, a MetaPath, where each path array are the different paths
  # that can be taken at each fork point
  def assemble_from(starting_path, visited_oriented_node_settables=Set.new)
    leash_length = @assembly_options[:max_bubble_length]
    current_bubble = nil
    if log.info? and starting_path.kind_of?(Bio::Velvet::Graph::OrientedNodeTrail)
      log.info "Assembling from: #{starting_path.to_shorthand}"
    end

    filterTips = lambda do |oneigh|
      is_tip, visiteds = is_short_tip?(oneigh)
      if is_tip
        visited_oriented_node_settables << oneigh.to_settable
        visiteds.each do |v|
          visited_oriented_node_settables << v
        end
      end
      log.debug "neighbour #{oneigh.to_shorthand} is_tip? #{is_tip}" if log.debug?
      is_tip
    end

    baseProblem = lambda do |oneigh|
      new_problem = DynamicProgrammingProblem.new
      new_problem.distance = 0
      new_path = Bio::Velvet::Graph::OrientedNodeTrail.new
      new_path.add_oriented_node oneigh
      new_problem.path = new_path
      new_problem.ubiquitous_oriented_nodes = Set.new
      new_problem.ubiquitous_oriented_nodes << oneigh.to_settable
      new_problem
    end

    extendedProblem = lambda do |problem, oneigh|
      new_problem = DynamicProgrammingProblem.new
      new_problem.distance = problem.distance + problem.path[-1].node.length_alone
      new_path = problem.path.copy
      new_path.add_oriented_node oneigh
      new_problem.path = new_path
      new_problem.ubiquitous_oriented_nodes = Set.new problem.ubiquitous_oriented_nodes
      new_problem.ubiquitous_oriented_nodes << oneigh.to_settable
      new_problem
    end

    metapath = MetaPath.new
    starting_path.each do |oriented_node|
      log.debug "adding onode at the start: #{oriented_node.to_shorthand}" if log.debug?
      metapath << oriented_node
    end

    # Keep track of nodes visited in this trajectory already so circuits can be avoided
    #visited_oriented_node_settables = Set.new
    starting_path.each do |e|
      if e.kind_of?(Bubble)
        e.oriented_nodes do |onode|
          visited_oriented_node_settables << onode.to_settable
        end
      else
        visited_oriented_node_settables << e.to_settable
      end
    end
    log.debug "Starting with visited nodes #{visited_oriented_node_settables.to_a.join(',')}" if log.debug?

    current_mode = :linear # :linear, :bubble, or :finished

    while current_mode != :finished
      if current_mode == :linear
        log.debug "Starting a non-bubble from #{metapath.to_shorthand}" if log.debug?
        while true
          oriented_neighbours = metapath.last_oriented_node.next_neighbours(@graph)
          log.debug "Found oriented neighbours #{oriented_neighbours.collect{|onode| onode.to_shorthand} }" if log.debug?

          legit_neighbours = nil
          # Cut off tips unless it is the only way
          if oriented_neighbours.length == 1
            legit_neighbours = oriented_neighbours
          else
            legit_neighbours = oriented_neighbours.reject &filterTips
          end

          if legit_neighbours.empty?
            # This is just a straight out dead end, and we can go no further.
            log.debug "Dead end reached" if log.debug?
            metapath.fate = MetaPath::DEAD_END_FATE
            current_mode = :finished
            break
          elsif legit_neighbours.length == 1
            # Linear thing here, just keep moving forward
            neighbour = legit_neighbours[0]

            # Stop if a circuit is detected
            # Tim - Always stop on a circuit in linear mode. Presumably means there is no way out.
            if visited_oriented_node_settables.include?(neighbour.to_settable)
              log.debug "Detected regular circuit by running into #{neighbour.to_settable}" if log.debug?
              metapath.fate = MetaPath::CIRCUIT_FATE
              current_mode = :finished
              break
            else
              visited_oriented_node_settables << neighbour.to_settable
              metapath << neighbour
            end

          else
            # Reached a fork in the graph here, the point of this algorithm, really.
            current_bubble = Bubble.new
            log.debug "Starting a bubble forking from metapath #{metapath.to_shorthand}" if log.debug?

            legit_neighbours.each do |oneigh|
              new_problem = baseProblem.call oneigh
              log.debug "Adding problem to bubble: #{new_problem}" if log.debug?

              current_bubble.enqueue new_problem
              current_mode = :bubble
            end
            break
          end
        end


      elsif current_mode == :bubble
        # We are in a bubble. Go get some.
        log.debug "entering bubble mode" if log.debug?

        # next problem = queue.shift. while distance of next problem is not beyond the leash length
        while current_mode == :bubble
          problem = current_bubble.shift

          if problem.nil?
            # Getting here seems improbable if not impossible.
            # The current bubble doesn't converge and just has short tips at the end, don't add it to the metapath
            metapath.fate = MetaPath::DEAD_END_FATE
            current_mode = :finished
            log.debug "Reached a dead end, ignoring this path" if log.debug?
            break
          end

          log.debug "Dequeued #{problem.to_shorthand}" if log.debug?
          if !leash_length.nil? and problem.distance > leash_length
            # The current bubble doesn't converge, don't add it to the metapath
            metapath.fate = MetaPath::DIVERGES_FATE
            current_mode = :finished
            log.debug "Bubble is past the leash length of #{leash_length}, giving up" if log.debug?
            break
          elsif current_bubble.convergent_on?(problem)
            log.debug "Bubble #{current_bubble.to_shorthand} convergent on #{problem.to_shorthand}" if log.debug?
            current_bubble.converge_on problem
            # convergement!
            # Bubble ended in a convergent fashion

            metapath << current_bubble
            # Add the nodes in the bubble to the list of visited nodes
            current_bubble.oriented_nodes do |onode|
              visited_oriented_node_settables << onode.to_settable
            end

            current_bubble = nil
            current_mode = :linear
            break
          else
            # otherwise we must search on in the bubble
            # get all neighbours that are not short tips
            log.debug "Bubble not convergent on #{problem.to_shorthand}" if log.debug?

            neighbours = problem.path.neighbours_of_last_node(@graph)

            # If there is only 1 way to go, go there
            if neighbours.length == 1
              log.debug "Only one way to go from this node, going there" if log.debug?

              oneigh = neighbours[0]
              new_problem = extendedProblem.call problem, oneigh
              new_problem.circular_path_detected = true if current_bubble.visited_oriented_nodes(problem).include? oneigh.to_settable
              current_bubble.enqueue new_problem
              log.debug "Enqueued #{new_problem.to_shorthand}, total nodes now #{current_bubble.num_known_problems} and num forks #{current_bubble.num_legit_forks}" if log.debug?

              # check to make sure we aren't going overboard in the bubbly-ness
              if !@assembly_options[:bubble_node_count_limit].nil? and current_bubble.num_known_problems > @assembly_options[:bubble_node_count_limit]
                log.debug "Too complex a bubble detected, giving up" if log.debug?
                metapath.fate = MetaPath::NODE_COUNT_LIMIT_REACHED
                current_mode = :finished
                break
              end
            else

              legit_neighbours = neighbours.reject &filterTips

              if legit_neighbours.length == 0
                # this is a kind of 'long' tip, possibly unlikely to happen much.
                # Forget about it and progress to the next problem having effectively
                # removed it from the bubble
                log.debug "Found no neighbours to re-enqueue" if log.debug?
              else
                # Increment complexity counter if this is a real fork
                if legit_neighbours.length > 1
                  current_bubble.num_legit_forks += 1
                end

                legit_neighbours.each do |oneigh|
                  new_problem = extendedProblem.call problem, oneigh
                  new_problem.circular_path_detected = true if current_bubble.visited_oriented_nodes(problem).include? oneigh.to_settable
                  current_bubble.enqueue new_problem
                  log.debug "Enqueued #{new_problem.to_shorthand}, total nodes now #{current_bubble.num_known_problems} and num forks #{current_bubble.num_legit_forks}" if log.debug?

                  # check to make sure we aren't going overboard in the bubbly-ness
                  if (!@assembly_options[:bubble_fork_limit].nil? and current_bubble.num_legit_forks > @assembly_options[:bubble_fork_limit]) or
                      (!@assembly_options[:bubble_node_count_limit].nil? and current_bubble.num_known_problems > @assembly_options[:bubble_node_count_limit])
                    log.debug "Too complex a bubble detected, giving up" if log.debug?
                    metapath.fate = MetaPath::NODE_COUNT_LIMIT_REACHED
                    current_mode = :finished
                    break
                  end
                end
              end
            end
          end
        end
      else
        raise "Programming error: Unexpected mode: #{current_mode}"
      end

      log.debug "Reached end of main loop in mode #{current_mode}" if log.debug?
    end

    return metapath, visited_oriented_node_settables
  end

  def seen_last_in_path?(path, seen_nodes)
    last = path[-1]
    if last.kind_of?(Bubble)
      return remove_seen_nodes_from_end_of_path(path, seen_nodes).length < path.length
    else
      return seen_nodes.include?(path[-1].to_settable)
    end
  end


  def remove_seen_nodes_from_end_of_path(path, seen_nodes)
    log.debug "Removing from the end of the path #{path.to_shorthand} any nodes in set of length #{seen_nodes.length}" if log.debug?

    node_seen = lambda do |oriented_node|
      seen_nodes.include?([oriented_node.node_id, Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST]) or
      seen_nodes.include?([oriented_node.node_id, Bio::Velvet::Graph::OrientedNodeTrail::END_IS_FIRST])
    end

    while !path.empty?
      last_node_or_bubble_index = path.length-1
      last_node_or_bubble = path[last_node_or_bubble_index]

      delete = false
      if last_node_or_bubble.kind_of?(Bubble)
        last_node_or_bubble.oriented_nodes do |onode|
          if node_seen.call(onode)
            delete = true
            break
          end
        end
      else
        delete = node_seen.call(last_node_or_bubble)
      end

      if delete
        path.delete_at last_node_or_bubble_index
      else
        # Last node is not previously seen, chop no further.
        break
      end
    end

    return path
  end

  def is_short_tip?(element)
    if element.kind_of?(Bubble)
      return false
    else
      super
    end
  end



  class MetaPath
    DIVERGES_FATE = 'diverges'
    DEAD_END_FATE = 'dead end'
    CIRCUIT_FATE = 'circuit'
    NODE_COUNT_LIMIT_REACHED = 'too many nodes in bubble'
    #CIRCUIT_WITHIN_BUBBLE_FATE = 'circuit within bubble' #Tim - shouldn't end metapath

    # How does this metapath end?
    attr_accessor :fate

    include Enumerable

    def initialize
      @internal_array = []
    end

    def each
      @internal_array.each do |e|
        yield e
      end
    end

    def [](index)
      @internal_array[index]
    end

    def delete_at(index)
      @internal_array.delete_at index
    end

    def empty?
      @internal_array.empty?
    end

    def last_oriented_node
      e = @internal_array[-1]
      if e.kind_of?(Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode)
        return e
      else
        # it is a bubble
        return e.converging_oriented_node
      end
    end

    def <<(oriented_node_or_bubble)
      @internal_array << oriented_node_or_bubble
    end
    alias_method :push, :<<

    def to_shorthand
      @internal_array.collect{|e| e.to_shorthand}.join(',')
    end

    def reverse!
      # Do regular reversal
      @internal_array.reverse!

      # Reverse all the internal parts
      @internal_array.each do |e|
        e.reverse!
      end

      return nil
    end

    def length
      @internal_array.length
    end

        # Yield all oriented nodes anywhere in the regular or bubble
    # bits.
    def each_oriented_node
      @internal_array.each do |e|
        if e.kind_of?(Bio::AssemblyGraphAlgorithms::BubblyAssembler::Bubble)
          e.oriented_nodes.each do |onode|
            yield onode
          end
        else
          yield e
        end
      end
    end

    def length_in_bp
      sum = 0
      each do |e|
        if e.kind_of?(Bio::AssemblyGraphAlgorithms::BubblyAssembler::Bubble)
          sum += e.reference_trail.length_in_bp_within_path
        else
          sum += e.node.length_alone
        end
      end
      return sum
    end

    def reference_trail
      trail = Bio::Velvet::Graph::OrientedNodeTrail.new

      trail.trail = collect do |e|
        if e.kind_of?(Bio::AssemblyGraphAlgorithms::BubblyAssembler::Bubble)
          e.reference_trail.trail
        else
          e
        end
      end.flatten

      return trail
    end

    def sequence
      reference_trail.sequence
    end

    def coverage
      coverages = []
      lengths = []
      each do |onode_or_bubble|
        if onode_or_bubble.kind_of?(Bio::AssemblyGraphAlgorithms::BubblyAssembler::Bubble)
          # Length isn't obvious, but let's go with reference path length just coz that's easy
          this_length = onode_or_bubble.reference_trail.length_in_bp_within_path
          lengths.push this_length

          # Coverage of a bubble is the coverage of each node in the bubble
          # each weighted by their length
          coverages.push onode_or_bubble.coverage
        else
          #regular node. So simple average coverage
          coverages.push onode_or_bubble.node.coverage
          lengths.push onode_or_bubble.node.length_alone
        end
      end

      # Then a simple weighted average
      i = -1
      total_length = lengths.reduce(:+)

      answer =  coverages.reduce(0.0) do |sum, cov|
        i += 1
        sum + (cov * lengths[i].to_f / total_length)
      end
      answer
    end
  end



  class Bubble
    include Bio::FinishM::Logging

    # The DynamicProgrammingProblem this bubble converges on
    attr_reader :converging_oriented_node_settable, :is_reverse

    # how many legit forks have been explored
    attr_accessor :num_legit_forks

    def initialize(options = {})
      @queue = DS::AnyPriorityQueue.new {|a,b| a<=b}
      @known_problems = {}
      @current_problems = Set.new
      @num_legit_forks = 0
      @max_cycles = options[:max_cycles] || DEFAULT_MAX_CYCLES
    end

    # Return the next closest dynamic programming problem,
    # removing it from the bubble
    def shift
      prob = @queue.shift
      @current_problems.delete prob.to_settable unless prob.nil?
      return prob
    end

    def visited_oriented_nodes(prob)
      #merge all visited nodes for relevant problems
      @known_problems[prob.to_settable].reduce(prob.ubiquitous_oriented_nodes) do |memo, problem|
        memo + problem.ubiquitous_oriented_nodes
      end
    end

    def ubiquitous_oriented_nodes(prob)
      #merge ubiquitous nodes for relevant problems
      @known_problems[prob.to_settable].reduce(prob.ubiquitous_oriented_nodes) do |memo, problem|
        memo & problem.ubiquitous_oriented_nodes
      end
    end

    def shortest_problem_distance(prob)
      # prioritise by the shortest distance for current problem
      @known_problems[prob.to_settable].collect{|prob| prob.distance}.min
    end


    def enqueue(dynamic_programming_problem)
      settable = dynamic_programming_problem.to_settable


      @known_problems[settable] ||= []
      @known_problems[settable].push dynamic_programming_problem

      # don't requeue current problem or circular problem
      unless dynamic_programming_problem.circular_path_detected == true or @current_problems.include? settable
        @queue.enqueue dynamic_programming_problem, shortest_problem_distance(dynamic_programming_problem)
        @current_problems << settable
      end
    end


    # return true if the given problem converges the bubble, else false
    def convergent_on?(dynamic_programming_problem)
      settable =  dynamic_programming_problem.to_settable

      @queue.each do |problem| #convergent until not
        return false unless ubiquitous_oriented_nodes(problem).include? settable
      end
      return true
    end

    # Finish off the bubble, assuming convergent_on? the given problem == true
    def converge_on(dynamic_programming_problem)
      @converging_oriented_node_settable = dynamic_programming_problem.to_settable
      #free some memory
      @queue = nil
      @current_problems = nil
    end

    # yield or failing that return an Array of the list of oriented_nodes found
    # in at least one path in this (presumed converged) bubble
    def oriented_nodes
      raise unless converged?
      seen_nodes = {}
      stack = DS::Stack.new
      initial_solution = @known_problems[@converging_oriented_node_settable][0]
      converging_onode = initial_solution.path[-1]
      stack.push converging_onode

      while onode = stack.pop
        settable = onode.to_settable
        next if seen_nodes.key?(settable)

        if block_given?
          if @is_reverse
            yield onode.reverse
          else
            yield onode
          end
        end

        seen_nodes[settable] = onode

        # queue neighbours for paths that don't contain the converging onode
        @known_problems[settable].each do |dpp|
          stack.push dpp.path[-2] unless dpp.path.length < 2 or dpp.path[0...-1].include? converging_onode
        end
      end

      return nil if block_given?
      return seen_nodes.values
    end

    def num_known_problems
      @known_problems.length
    end

    # Iterate over the paths returning each as an OrientedNodeTrail.
    # Assumes the path is convergent.
    def each_path(max_cycles = @max_cycles)
      raise unless converged?

      log.debug "Iterating through each path of bubble" if log.debug?

      # Tim - use priority queue to yield shortest paths first
      queue = DS::AnyPriorityQueue.new {|a, b| a<=b}
      counter = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder::CycleCounter.new max_cycles
      initial_solution = @known_problems[@converging_oriented_node_settable][0]
      queue.enqueue [initial_solution.path, [], 0], 0
      #log.debug "Pushed to stack #{initial_solution.path.to_shorthand}" if log.debug?

      while path_parts = queue.dequeue
        direct_node_trail = path_parts[0]
        second_part = path_parts[1]
        second_distance = path_parts[2]
        #log.debug "Popped #{direct_node_trail.to_shorthand} and [#{second_part.collect{|o| o.to_shorthand}.join(',') }]" if log.debug?


        if direct_node_trail.trail.length == 0
          yield_path = Bio::Velvet::Graph::OrientedNodeTrail.new
          yield_path.trail = second_part
          if @is_reverse
            yield yield_path.reverse
          else
            #log.debug "Yielded #{yield_path.to_shorthand}" if log.debug?
            yield yield_path
          end
        else
          # go down the path, looking for other paths
          head_onode = direct_node_trail.trail[-1]
          new_second_part = [head_onode]+second_part
          new_second_distance = second_distance+head_onode.node.length_alone
          if second_part.length > 1 and head_onode.to_settable == @converging_oriented_node_settable
            log.debug "Ignoring path with cycle through converged node." if log.debug?
            next
          end
          if second_part.include? head_onode
            log.debug "Cycle at node #{head_onode.node_id} in path #{second_part.collect{|onode| onode.node.node_id}.join(',')}." if log.debug?
            @circuitous = true unless @circuitous
            if max_cycles == 0 or max_cycles < counter.path_cycle_count(new_second_part)
              log.debug "Not finishing cyclic path with too many cycles." if log.debug?
              next
            end
          end

          new_problems = @known_problems[head_onode.to_settable]
          #log.debug "Found new problems: #{new_problems.collect{|prob| prob.to_shorthand}.join(' ') }" if log.debug?
          problem_leads = Set.new
          new_problems.each do |new_problem|
            # Only enqueue paths where the second-to-head onode is not already queued
            unless new_problem.path.length < 2
              lead_settable = new_problem.path[-2].to_settable
              next if problem_leads.include? lead_settable
              problem_leads << lead_settable
            end

            # TODO: deal with circuits
            new_trail = Bio::Velvet::Graph::OrientedNodeTrail.new
            new_trail.trail = new_problem.path[0...-1]
            #log.debug "Enqueuing #{new_trail.to_shorthand} and [#{new_second_part.collect{|o| o.to_shorthand}.join(',') }]" if log.debug?
            queue.enqueue [new_trail, new_second_part, new_second_distance], new_second_distance
          end
        end
      end
    end

    def paths
      to_return = []
      each_path do |path|
        to_return.push path
      end
      to_return
    end

    def converged?
      !@converging_oriented_node_settable.nil?
    end

    # Return the OrientedNode that converges this bubble, behaviour
    # undefined if bubble is not converged
    def converging_oriented_node
      @known_problems[@converging_oriented_node_settable][0].path[-1]
    end

    def to_shorthand
      shorts = []
      if converged?
        shorts = paths.sort{|a,b| a.to_shorthand <=> b.to_shorthand }.collect{|path| path.to_shorthand}
      else
        @queue.each do |problem|
          shorts.push problem.to_shorthand
        end
      end
      return "{#{shorts.join('|') }}"
    end

    def reverse!
      @is_reverse ||= false
      @is_reverse = !@is_reverse
    end

    # This doesn't make sense unless this is a converged bubble and the index == -1
    # because otherwise there is multiple answers
    def [](index)
      raise unless index == -1
      return Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(
        @converging_oriented_node_settable[0],
        @converging_oriented_node_settable[1]
        )
    end

    # Return one trail that exemplifies the paths through this bubble.
    # Current method of path selection is simply greedy, taking the highest coverage node
    # at each fork (or failing that the node with the lower node_id).
    def reference_trail(max_cycles = @max_cycles)
      raise unless converged?

      log.debug "Finding reference trail" if log.debug?

      reference_trail = []
      #counter = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder::CycleCounter.new max_cycles
      comparator = lambda do |problem1, problem2|
        onode1 = nil
        onode2 = nil
        if problem1.path.length == 1 and problem2.path.length > 1
          # Here the comparison cannot be made on 2nd last node coverages
          # since one of the paths goes straight from the initial to the terminal
          # node. Choose instead based on if the second last node has higher or lower
          # coverage than the final node
          onode1 = problem1.path[-1]
          onode2 = problem2.path[-2]
        elsif problem2.path.length == 1 and problem1.path.length > 1
          onode1 = problem1.path[-2]
          onode2 = problem2.path[-1]
        else
          onode1 = problem1.path[-2]
          onode2 = problem2.path[-2]
        end
        log.debug "Comparing nodes #{onode1.node.node_id} and #{onode2.node.node_id}" if log.debug?

        # Prefer non-cyclic routes
        log.debug "Looking for nodes in reference trail #{reference_trail.collect{|onode| onode.node.node_id}.join(',')}" if log.debug?
        return -1 if reference_trail.include? onode1 and
            (problem2.path.length == 1  or !reference_trail.include? onode2)
        return 1 if reference_trail.include? onode2 and
            (problem1.path.length == 1 or !reference_trail.include? onode1)

        log.debug "Comparing coverages" if log.debug?
        if onode1.node.coverage == onode2.node.coverage
          log.debug "Coverage is the same" if log.debug?
          -(onode1.node.node_id <=> onode2.node.node_id)
        else
          onode1.node.coverage <=> onode2.node.coverage
        end
      end


      current_problem = @known_problems[@converging_oriented_node_settable].max do |problem1, problem2|
        comparator.call problem1, problem2
      end

      while !current_problem.nil?
        onode_to_add = current_problem.path[-1]
        log.debug "Preferred path #{current_problem.path.collect{|onode| onode.node.node_id}.join(',')}" if log.debug?
        #raise CircuitousPathDetected if reference_trail.include?(onode_to_add)
        reference_trail.push onode_to_add

        if current_problem.path.length == 1
          current_problem = nil
        else
          current_problem = @known_problems[current_problem.path[-2].to_settable].max do |problem1, problem2|
            comparator.call problem1, problem2
          end
        end
      end

      node_trail = Bio::Velvet::Graph::OrientedNodeTrail.new
      node_trail.trail = reference_trail.reverse
      node_trail = node_trail.reverse if @is_reverse

      return node_trail
    end

    # Does this (coverged) bubble contain any circuits?
    def circuitous?
      raise unless converged?
      each_path(0) {|| break if @circuitous}
      @circuitous == true
    end

    # Coverage of a bubble is the coverage of each node in the bubble
    # each weighted by their length
    def coverage
      sum = 0.0
      length = 0
      oriented_nodes do |onode|
        node_length = onode.node.length_alone
        sum += onode.node.coverage * node_length
        length += node_length
      end
      return sum / length
    end
  end

  class DynamicProgrammingProblem
    attr_accessor :path, :ubiquitous_oriented_nodes, :distance, :circular_path_detected

    def initialize
      @path = []
      @ubiquitous_oriented_nodes = Set.new
    end

    def to_settable
      @path[-1].to_settable
    end

    def to_s
      ubiquitous_nodes = @ubiquitous_oriented_nodes.collect do |settabled|
        "#{settabled[0] }#{settabled[1] == Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST ? 's' : 'e'}"
      end
      return "DPP #{self.object_id}: #{@path.to_shorthand}/#{ubiquitous_nodes.join(',') }/#{distance}"
    end

    def to_shorthand
      ubiquitous_nodes = @ubiquitous_oriented_nodes.collect do |settabled|
        "#{settabled[0] }#{settabled[1] == Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST ? 's' : 'e'}"
      end
      "#{@path.to_shorthand}/#{ubiquitous_nodes.join(',') }/#{distance}"
    end
  end

  class CircuitousPathDetected < Exception; end
end
