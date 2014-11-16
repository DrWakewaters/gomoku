/*
	if(isBlack && selfHasAWinningPly) {
		board->setBlackAtPosition(winningPositionSelf);
		if(!AI::gameWon(board, board->black)) {
			exit(-1);
		}
		board->setEmptyAtPosition(winningPositionSelf);
	}
	if(isBlack && opponentHasAWinningPly) {		
		board->setWhiteAtPosition(winningPositionOpponent);
		if(!AI::gameWon(board, board->white)) {
			exit(-1);
		}
		board->setEmptyAtPosition(winningPositionOpponent);
	}
	if(!isBlack && selfHasAWinningPly) {		
		board->setWhiteAtPosition(winningPositionSelf);
		if(!AI::gameWon(board, board->white)) {
			exit(-1);
		}
		board->setEmptyAtPosition(winningPositionSelf);
	}
	if(!isBlack && opponentHasAWinningPly) {		
		board->setBlackAtPosition(winningPositionOpponent);
		if(!AI::gameWon(board, board->black)) {
			exit(-1);
		}
		board->setEmptyAtPosition(winningPositionOpponent);
	}
*/

#include "AI.h"

int computed = 0;
int computedForMoveOrder = 0;
int hashes = 0;

// Compute the score using iterative deepening and the alpha-beta pruning algorithm.
Position AI::computeScore(Board *board, int8_t maxSearchDepth, double maxSearchTime, bool isBlack, bool printAIInformation) {
	// 
	std::vector<std::vector<Position>> nextPositions = std::vector<std::vector<Position>>(maxSearchDepth, std::vector<Position>(board->boardWidth*board->boardHeight));
	// This is the iterative deepening loop. Its purpose is to speed up the alphabeta algorithm. It does this by finding out a good order in which to try plies, which makes for more aggressive pruning.
	auto startTime = std::chrono::high_resolution_clock::now();
	double computationTime = 0.0;
	// This is the iterative deepening loop. Its purpose is to speed up the alphabeta algorithm. It does this by finding the best next ply for most board positions. At the next iteration of the deepener, this promising ply (it is not certain that it is the best ply now that the search depth is one level deeper, but it often is) is tried first. This gives aggressive pruning.
	std::vector<Position> principalVariation;
	int16_t score;
	bool gameWon = false;
	for(int8_t i = 1; i <= maxSearchDepth && computationTime < maxSearchTime; i++) {
		AI::computeScoreUsingAlphaBeta(board, nextPositions, i, i, MIN_SCORE, MAX_SCORE, isBlack, false);
		principalVariation = board->principalVariation(isBlack, i, &score, &gameWon);
		if(printAIInformation) {
			std::cout << "The evaluation method was called " << computed << " times (for move order: " << computedForMoveOrder << ", at leaves: " << computed-computedForMoveOrder << ")."  << std::endl; 
			std::cout << "The hash table has " << hashes << " hash objects of size " << sizeof(Hash) << " bytes each." << std::endl;
			if(principalVariation.size() > 0) {
				std::cout << "The principal variation: ";
				for(size_t j = 0; j < principalVariation.size()-1; j++) {
					std::cout << "(" << static_cast<int>(principalVariation.at(j).y) << ", " << static_cast<int>(principalVariation.at(j).x) << "), ";
				}
				std::cout << "(" << static_cast<int>(principalVariation.back().y) << ", " << static_cast<int>(principalVariation.back().x) << ")." << std::endl;
				std::cout << "The score afterwards is " << static_cast<int>(score) << "." << std::endl;
			}
			std::cout << std::endl;
		}
		board->clearHashTablePartially();
		computed = 0;
		computedForMoveOrder = 0;
		auto endTime = std::chrono::high_resolution_clock::now();
		computationTime = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(endTime-startTime).count())/1000000000.0;
		// If we know that we can win or if we know that the opponent can win if it searches at least to this depth (of course, the latter might not be the case, but I do not think the algorithm will switch to another ply if we search deeper), do not search deeper.
		if(gameWon) {
			break;
		}
	}
	hashes = 0;
	return principalVariation.at(0);
}

int16_t AI::computeScoreUsingAlphaBeta(Board *board, std::vector<std::vector<Position>> &nextPositions, int8_t maxLevels, int8_t levelsToBottom, int16_t alpha, int16_t beta, bool isBlack, bool mustComputeScore) {	
	// The variable hashIndexPrimary is the hash value for this position, modulo the size of the hash table. There might be several different boards with the same hashIndexPrimary, therefore a second index, hashIndexSecondary, is used.
	int hashIndexPrimary = static_cast<int>((board->hashValue)%(board->hashTable.size()));
	// Either return the inner index of an already existing hash object (that corresponds to this board position), or create a new one.
	int hashIndexSecondary = AI::setupHashObject(board, hashIndexPrimary);
	// If we have a valid score for this hash, then return it.
	if(board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).hasValidScore) {
		return board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).score;
	}
	// If either side has won, the score is maximal or minimal.
	if(isBlack && AI::gameWon(board, board->white)) {
		board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).hasValidScore = true;
		board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).score = MIN_SCORE;
		return MIN_SCORE;
	} else if(!isBlack && AI::gameWon(board, board->black)) {
		board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).hasValidScore = true;
		board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).score = MAX_SCORE;
		return MAX_SCORE;
	}
	// If we are at a leaf position (that is, we will not search deeper), evaluate the position and return its score.
	if(levelsToBottom == 0 || board->boardIsFull() || mustComputeScore) {
		int16_t score = AI::evaluateBoard(board, isBlack, false);		
		board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).hasValidScore = true;
		board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).score = score;
		return score;
	}
	// If we have a half open (or open) four, we should make the winning ply. If not and if the opponent has a half open (or open) four, we must stop the attack.
	// The same goes if the we have or the opponent has a 3-empty-1 or a 2-empty-2.
	// In either case, when alphaBetaInnerLoop is called, it will be with just this ply (if there is one ply for both us and the opponent, our ply will be chosen).
	// If the opponent has a half open four AND we do not have a half open four AND we are on the top level, then we tell computeScoreUsingAlphaBeta to compute the score the next time it is called.
	Position winningPositionSelf = Position();
	Position winningPositionOpponent = Position();
	bool selfHasAWinningPly = false;
	bool opponentHasAWinningPly = false;
	if(isBlack) {
		selfHasAWinningPly = hasAWinningPly(board, board->black, board->empty, winningPositionSelf);
		opponentHasAWinningPly = hasAWinningPly(board, board->white, board->empty, winningPositionOpponent);
	} else {
		selfHasAWinningPly = hasAWinningPly(board, board->white, board->empty, winningPositionSelf);
		opponentHasAWinningPly = hasAWinningPly(board, board->black, board->empty, winningPositionOpponent);
	}
	if(opponentHasAWinningPly && !selfHasAWinningPly && (levelsToBottom == maxLevels)) {
		mustComputeScore = true;
	}
	// Find the best ply. Then store its position in the hash object (to save memory, the score is not saved).
	Ply bestNextPly;
	if(isBlack) {
		bestNextPly.score = MIN_SCORE;
	} else {
		bestNextPly.score = MAX_SCORE;		
	}
	bool hasFoundABestNextPly = false;	
	// We now begin the main loop in the alpha beta algorithm. We check plies, which will modify alpha (if we are the maximizing player) or beta (if we are the minimizing player).
	// If this causes beta <= alpha, then we do not check any more plies, since they are not relevant (we prune that remaining part of the search tree).
	// We first check the most promising ply, then (if that doesn't force an alpha-beta pruning, which we hope it does) the other possible plies (until we've checked them all or until there is a pruning). How did we find the promising ply? At earlier turns in the iterative deepener (ie with maxLevels smaller than it is now) it was found to be the best one. This information was stored and is now used (and at this level we store the same kind of information for the next turn in the iterative deepener).
	if(selfHasAWinningPly) {
		AI::alphaBetaInnerLoop(board, nextPositions, winningPositionSelf, &alpha, &beta, levelsToBottom, maxLevels, isBlack, mustComputeScore, &bestNextPly, &hasFoundABestNextPly);
	} else if(opponentHasAWinningPly) {
		AI::alphaBetaInnerLoop(board, nextPositions, winningPositionOpponent, &alpha, &beta, levelsToBottom, maxLevels, isBlack, maxLevels == levelsToBottom, &bestNextPly, &hasFoundABestNextPly);		
	} else {
		Position oldBestNextPosition = board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).bestNextPosition;
		if(board->interestingAtPosition(oldBestNextPosition)) {
			AI::alphaBetaInnerLoop(board, nextPositions, oldBestNextPosition, &alpha, &beta, levelsToBottom, maxLevels, isBlack, mustComputeScore, &bestNextPly, &hasFoundABestNextPly);
			if(beta <= alpha) {
				goto alphabetaprune;
			}
		}
		AI::computeNextPositions(board, nextPositions, hashIndexPrimary, hashIndexSecondary, maxLevels, levelsToBottom, isBlack);
		for(Position nextPosition : nextPositions.at(maxLevels-levelsToBottom)) {
			AI::alphaBetaInnerLoop(board, nextPositions, nextPosition, &alpha, &beta, levelsToBottom, maxLevels, isBlack, mustComputeScore, &bestNextPly, &hasFoundABestNextPly);
			if(beta <= alpha) {
				goto alphabetaprune;
			}
		}
	}
	alphabetaprune:
	// Store the best ply. (There could be a position stored from an earlier iteration in the deepener that is replaced - this is generally good. There could also be a position from this iteration that is replaced - this can be good or bad).
	board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).bestNextPositionLevel = maxLevels;
	board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).bestNextPosition = bestNextPly.position;
	// If we are the maximizing player, we return alpha, otherwise beta. Hey, remember to use cutoff.
	if(isBlack) {
		return alpha;
	}
	return beta;

}

void AI::computeNextPositions(Board *board, std::vector<std::vector<Position>> &nextPositions, int hashIndexPrimary, int hashIndexSecondary, int8_t maxLevels, int8_t levelsToBottom, bool isBlack) {
	Position oldBestNextPosition = board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).bestNextPosition;
	nextPositions.at(maxLevels-levelsToBottom).resize(0);
	Rectangle searchRectangle = board->searchRectangle;
	for(int i = searchRectangle.upperLeft.y; i <= searchRectangle.lowerRight.y; i++) {
		for(int j = searchRectangle.upperLeft.x; j <= searchRectangle.lowerRight.x; j++) {
			Position position = Position(i, j);
			if(!board->interestingAtPosition(position)) {
				continue;
			}
			if(position.y == oldBestNextPosition.y && position.x == oldBestNextPosition.y) {
				continue;
			}
			nextPositions.at(maxLevels-levelsToBottom).push_back(position);
		}
	}
	// Unless we are one the level above the leaves, compute the scores for each possible next ply (by making it, evaluating it and then unmaking it) and then sort them (keeping the old best next position, computed earlier in the iterative deepener, first!).
	if(levelsToBottom < 2) {
		return;
	}
	size_t numberOfNextPositions = nextPositions.at(maxLevels-levelsToBottom).size();
	std::vector<int> scores(numberOfNextPositions);
	for(size_t i = 1; i < numberOfNextPositions; i++) {
		// Store the interesting-vector and the search rectangle; they will be restored later.	
		std::vector<int> interestingBackup(board->interesting);
		Rectangle searchRectangleBackup(board->searchRectangle);
		// Make the ply and modify the interesting-vector.
		if(isBlack) {
			board->setBlackAtPosition(nextPositions.at(maxLevels-levelsToBottom).at(i));
		} else {
			board->setWhiteAtPosition(nextPositions.at(maxLevels-levelsToBottom).at(i));
		}
		board->updateInteresting(nextPositions.at(maxLevels-levelsToBottom).at(i));
		// Compute the score.
		scores.at(i) = AI::evaluateBoard(board, isBlack, false);
		computedForMoveOrder++;
		// Undo the ply and restore the vector with interesting plies and the search rectangle.
		board->setEmptyAtPosition(nextPositions.at(maxLevels-levelsToBottom).at(i));
		board->interesting = interestingBackup;
		board->searchRectangle = searchRectangleBackup;
	}
	// If we are black, we should sort the plies from highest score to lowest; otherwise lowest to highest. If start == 1, then there is a position from the iterative deepener in the beginning: do not move it.
	if(isBlack) {
		for(size_t i = 0; i < numberOfNextPositions-1; i++) {
			for(size_t j = i+1; j > 0; j--) {
				if(scores.at(j-1) < scores.at(j)) {
					int tmpScore = scores.at(j-1);
					scores.at(j-1) = scores.at(j);
					scores.at(j) = tmpScore;
					Position tmpPosition = nextPositions.at(maxLevels-levelsToBottom).at(j-1);
					nextPositions.at(maxLevels-levelsToBottom).at(j-1) = nextPositions.at(maxLevels-levelsToBottom).at(j);
					nextPositions.at(maxLevels-levelsToBottom).at(j) = tmpPosition;
				}
			}
		}
	} else {
		for(size_t i = 0; i < numberOfNextPositions-1; i++) {
			for(size_t j = i+1; j > 0; j--) {
				if(scores.at(j-1) > scores.at(j)) {
					int tmpScore = scores.at(j-1);
					scores.at(j-1) = scores.at(j);
					scores.at(j) = tmpScore;
					Position tmpPosition = nextPositions.at(maxLevels-levelsToBottom).at(j-1);
					nextPositions.at(maxLevels-levelsToBottom).at(j-1) = nextPositions.at(maxLevels-levelsToBottom).at(j);
					nextPositions.at(maxLevels-levelsToBottom).at(j) = tmpPosition;
				}
			}
		}
	}
}

int AI::setupHashObject(Board *board, int hashIndexPrimary) {
	// This is the 3*16 rightmost bits of the hash value of this board position (only those bits are stored in a hash object).
	uint16_t hashValue[3] = {static_cast<uint16_t>(board->hashValue >> 48), static_cast<uint16_t>((board->hashValue >> 32)%(1 << 16)), static_cast<uint16_t>((board->hashValue >> 16)%(1 << 16))};
	// First check if we have a hash object for this board position.
	for(int i = 0; i < static_cast<int>(board->hashTable.at(hashIndexPrimary).size()); i++) {
		if(board->hashTable.at(hashIndexPrimary).at(i).hashValue[0] == hashValue[0] && board->hashTable.at(hashIndexPrimary).at(i).hashValue[1] == hashValue[1] && board->hashTable.at(hashIndexPrimary).at(i).hashValue[2] == hashValue[2]) {
			return i;
		}
	}
	// Ok, so we do not have a hash object for this position. Therefore we create one.
	board->hashTable.at(hashIndexPrimary).push_back(Hash());
	hashes++;
	int hashIndexSecondary = static_cast<int>(board->hashTable.at(hashIndexPrimary).size()) - 1;
	board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).hashValue[0] = hashValue[0];
	board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).hashValue[1] = hashValue[1];
	board->hashTable.at(hashIndexPrimary).at(hashIndexSecondary).hashValue[2] = hashValue[2];
	return hashIndexSecondary;
}

void AI::alphaBetaInnerLoop(Board *board, std::vector<std::vector<Position>> &nextPositions, Position nextPosition, int16_t *alpha, int16_t *beta, int8_t levelsToBottom, int8_t maxLevels, bool isBlack, bool mustComputeScore, Ply *bestNextPly, bool *hasFoundABestNextPly) {
	// Store the interesting-vector and the search rectangle; they will be restored later.	
	std::vector<int> interestingBackup(board->interesting);
	Rectangle searchRectangleBackup(board->searchRectangle);
	// Make the ply and modify the interesting-vector.
	if(isBlack) {
		board->setBlackAtPosition(nextPosition);
	} else {
		board->setWhiteAtPosition(nextPosition);
	}
	board->updateInteresting(nextPosition);
	// Compute the score recursively by calling computeScoreUsingAlphaBeta for one level deeper.
	int16_t score = AI::computeScoreUsingAlphaBeta(board, nextPositions, maxLevels, levelsToBottom-1, *alpha, *beta, !isBlack, mustComputeScore);
	if((isBlack && score > bestNextPly->score) || (!isBlack && score < bestNextPly->score) || !*hasFoundABestNextPly) {
		*hasFoundABestNextPly = true;
		bestNextPly->score = score;
		bestNextPly->position = nextPosition;
	}
	// Undo the ply and restore the vector with interesting plies and the search rectangle.
	board->setEmptyAtPosition(nextPosition);
	board->interesting = interestingBackup;
	board->searchRectangle = searchRectangleBackup;
	// If we have found a new best value, store it it either alpha or beta. If we also happen are at the top of the search tree (at the first ply, the one that we will actually do), store the ply.
	if(isBlack && score > *alpha) {
		*alpha = score;
	} else if(!isBlack && score < *beta) {
		*beta = score;
	}
}

// Evaluate the board by counting the number of n-in-a-row:s the players have. If an n-in-a-row is unblocked, it is worth more than if it is blocked on one side. It is worth zero if it is blocked on both sides. An n-in-a-row is worth much more if we are the person who will do the nex ply (since we are then guaranteed to immediately get an n+1-in-a-row, if we wish). This evaluation method is very simple, and probably not that good - but it works ok.
int16_t AI::evaluateBoard(Board *board, bool isBlack, bool printInformation) {
	computed++;
	double fivesBlack = 0.0;
	double foursBlackOpen = static_cast<double>(AI::numberOfFoursOpen(board, board->black, board->empty));
	double foursBlackHalfOpen = static_cast<double>(AI::numberOfFoursHalfOpen(board, board->black, board->empty)) - foursBlackOpen;
	double threesBlackOpen = static_cast<double>(AI::numberOfThreesOpen(board, board->black, board->empty));
	double threesBlackHalfOpen = static_cast<double>(AI::numberOfThreesHalfOpen(board, board->black, board->empty)) - 2*foursBlackOpen - foursBlackHalfOpen - threesBlackOpen;
	double twosBlackOpen = static_cast<double>(AI::numberOfTwosOpen(board, board->black, board->empty));
	double twosBlackHalfOpen = static_cast<double>(AI::numberOfTwosHalfOpen(board, board->black, board->empty)) - 2*foursBlackOpen - foursBlackHalfOpen - 2*threesBlackOpen - threesBlackHalfOpen - twosBlackOpen;
	double fivesWhite = 0.0;
	double foursWhiteOpen = static_cast<double>(AI::numberOfFoursOpen(board, board->white, board->empty));
	double foursWhiteHalfOpen = static_cast<double>(AI::numberOfFoursHalfOpen(board, board->white, board->empty)) - foursWhiteOpen;
	double threesWhiteOpen = static_cast<double>(AI::numberOfThreesOpen(board, board->white, board->empty));
	double threesWhiteHalfOpen = static_cast<double>(AI::numberOfThreesHalfOpen(board, board->white, board->empty)) - 2*foursWhiteOpen - foursWhiteHalfOpen - threesWhiteOpen;
	double twosWhiteOpen = static_cast<double>(AI::numberOfTwosOpen(board, board->white, board->empty));
	double twosWhiteHalfOpen = static_cast<double>(AI::numberOfTwosHalfOpen(board, board->white, board->empty)) - 2*foursWhiteOpen - foursWhiteHalfOpen - 2*threesWhiteOpen - threesWhiteHalfOpen - twosWhiteOpen;
	// The player who will make the next ply has an advantage. These modifiers tries to take that into account, but they are quite crude.
	if(isBlack) {
		if(foursBlackOpen > 0.0 || foursBlackHalfOpen > 0.0) {
			fivesBlack += 0.5;
		} else if(threesBlackOpen > 0.0) {
			foursBlackOpen += 0.5;
		} else if(threesBlackHalfOpen > 0.0) {
			foursBlackHalfOpen += 0.5;
		} else if(twosBlackOpen > 0.0) {
			threesBlackOpen += 0.5;
		} else if(twosBlackHalfOpen > 0.0) {
			threesBlackHalfOpen += 0.5;
		}
	} else {
		if(foursWhiteOpen > 0.0 || foursWhiteHalfOpen > 0.0) {
			fivesWhite += 0.5;
		} else if(threesWhiteOpen > 0.0) {
			foursWhiteOpen += 0.5;
		} else if(threesWhiteHalfOpen > 0.0) {
			foursWhiteHalfOpen += 0.5;
		} else if(twosWhiteOpen > 0.0) {
			threesWhiteOpen += 0.5;
		} else if(twosWhiteHalfOpen > 0.0) {
			threesWhiteHalfOpen += 0.5;
		}
	}
	if(printInformation) {
		std::cout << "Status right now (half a point means a possibility to get it after this ply)." << std::endl;
		std::cout << "Black, open (4 3 2): " << foursBlackOpen << " " << threesBlackOpen << " " << twosBlackOpen << "." << std::endl;
		std::cout << "Black, half open (4 3 2): " << foursBlackHalfOpen << " " << threesBlackHalfOpen << " " << twosBlackHalfOpen << std::endl;
		std::cout << "White, open (4 3 2): " << foursWhiteOpen << " " << threesWhiteOpen << " " << twosWhiteOpen << "." << std::endl;
		std::cout << "White, half open (4 3 2): " << foursWhiteHalfOpen << " " << threesWhiteHalfOpen << " " << twosWhiteHalfOpen << std::endl << std::endl;
	}
	// These numbers are quite random. A better evaluation method would give a much better AI, but is not trivial to implement. Perhaps AIs with different constants could battle each other (evolutionary programming perhaps?)?
	int blackScore = 4000.0*fivesBlack + 400.0*foursBlackOpen + 40.0*foursBlackHalfOpen + 20.0*threesBlackOpen + 10.0*threesBlackHalfOpen + 4.0*twosBlackOpen + 2.0*twosBlackHalfOpen;
	int whiteScore = 4000.0*fivesWhite + 400.0*foursWhiteOpen + 40.0*foursWhiteHalfOpen + 20.0*threesWhiteOpen + 10.0*threesWhiteHalfOpen + 4.0*twosWhiteOpen + 2.0*twosWhiteHalfOpen;
	return blackScore - whiteScore;
}

bool AI::gameWon(Board *board, bool isBlack) {
	if(isBlack) {
		return AI::gameWon(board, board->black);
	}
	return AI::gameWon(board, board->white);
}

bool AI::gameWon(Board *board, std::vector<int> &bitBoard) {
	int domain;
	int uppermostStone = board->searchRectangle.upperLeft.y;
	int lowermostStone = board->searchRectangle.lowerRight.y;
	for(int i = uppermostStone; i <= lowermostStone-4; i++) {
		// Vertically.
		domain = bitBoard[i] & bitBoard[i+1] & bitBoard[i+2] & bitBoard[i+3] & bitBoard[i+4];
		if(domain != 0) {
			return true;
		}
		// Diagonally \ (basically shift and then check vertically).
		domain = bitBoard[i] & (bitBoard[i+1] >> 1) & (bitBoard[i+2] >> 2) & (bitBoard[i+3] >> 3) & (bitBoard[i+4] >> 4);
		if(domain != 0) {
			return true;
		}
		// Diagonally / (basically shift and then check vertically).
		domain = bitBoard[i] & (bitBoard[i+1] << 1) & (bitBoard[i+2] << 2) & (bitBoard[i+3] << 3) & (bitBoard[i+4] << 4);
		if(domain != 0) {
			return true;
		}
	}
	// Horizontally.
	for(int i = uppermostStone; i <= lowermostStone; i++) {
		domain = bitBoard[i] & (bitBoard[i] >> 1) & (bitBoard[i] >> 2) & (bitBoard[i] >> 3) & (bitBoard[i] >> 4);
		if(domain != 0) {
			return true;
		}
	}
	return false;
}

bool AI::setwinningPosition(int leftColumnLimit, int rightColumnLimit, int row, int columnNumber, int columnModifier, Position &winningPosition) {
	if(columnNumber != 0) {
		for(int column = leftColumnLimit; column <= rightColumnLimit; column++) {
			if((columnNumber & (1 << column)) != 0) {
				winningPosition.y = static_cast<uint8_t>(row);
				winningPosition.x = static_cast<uint8_t>(column + columnModifier);
			}
		}
		return true;
	}
	return false;
}

bool AI::hasAWinningPly(Board *board, std::vector<int> &bitBoard, std::vector<int> &bitBoardEmpty, Position &winningPosition) {
	int domain;
	int uppermostStone = board->searchRectangle.upperLeft.y;
	int lowermostStone = board->searchRectangle.lowerRight.y;
	int leftColumnLimit = board->searchRectangle.upperLeft.x;
	int rightColumnLimit = board->searchRectangle.lowerRight.x;
	for(int i = uppermostStone+1; i <= lowermostStone-4; i++) {
		// Vertically.
		domain = (bitBoardEmpty[i-1] | bitBoardEmpty[i+4]) & bitBoard[i] & bitBoard[i+1] & bitBoard[i+2] & bitBoard[i+3];
		if(domain != 0) {
			if(!AI::setwinningPosition(leftColumnLimit, rightColumnLimit, i-1, domain & bitBoardEmpty[i-1], 0, winningPosition)) {
				AI::setwinningPosition(leftColumnLimit, rightColumnLimit, i+4, domain & bitBoardEmpty[i+4], 0, winningPosition);				
			}
			return true;
		}
		// Diagonally \ (basically shift and then check vertically).
		domain = ((bitBoardEmpty[i-1] << 1) | (bitBoardEmpty[i+4] >> 4)) & bitBoard[i] & (bitBoard[i+1] >> 1) & (bitBoard[i+2] >> 2) & (bitBoard[i+3] >> 3);
		if(domain != 0) {
			if(!AI::setwinningPosition(leftColumnLimit, rightColumnLimit, i-1, domain & (bitBoardEmpty[i-1] << 1), -1, winningPosition)) {
				AI::setwinningPosition(leftColumnLimit, rightColumnLimit, i+4, domain & (bitBoardEmpty[i+4] >> 4), 4, winningPosition);				
			}
			return true;
		}
		// Diagonally / (basically shift and then check vertically).
		domain = ((bitBoardEmpty[i-1] >> 1) | (bitBoardEmpty[i+4] << 4)) & bitBoard[i] & (bitBoard[i+1] << 1) & (bitBoard[i+2] << 2) & (bitBoard[i+3] << 3);
		if(domain != 0) {
			if(!AI::setwinningPosition(leftColumnLimit, rightColumnLimit, i-1, domain & (bitBoardEmpty[i-1] >> 1), 1, winningPosition)) {
				AI::setwinningPosition(leftColumnLimit, rightColumnLimit, i+4, domain & (bitBoardEmpty[i+4] << 4), -4, winningPosition);				
			}
			return true;
		}
	}
	// We have to treat the top and bottom separately.
	if(uppermostStone+4 <= lowermostStone) {
		// Vertically at the top.
		domain = bitBoard[uppermostStone] & bitBoard[uppermostStone+1] & bitBoard[uppermostStone+2] & bitBoard[uppermostStone+3] & bitBoardEmpty[uppermostStone+4];
		if(domain != 0) {
			AI::setwinningPosition(leftColumnLimit, rightColumnLimit, uppermostStone+4, domain & bitBoardEmpty[uppermostStone+4], 0, winningPosition);
			return true;
		}
		// Diagonally \ at the top (basically shift and then check vertically).
		domain = bitBoard[uppermostStone] & (bitBoard[uppermostStone+1] >> 1) & (bitBoard[uppermostStone+2] >> 2) & (bitBoard[uppermostStone+3] >> 3) & (bitBoardEmpty[uppermostStone+4] >> 4);
		if(domain != 0) {
			AI::setwinningPosition(leftColumnLimit, rightColumnLimit, uppermostStone+4, domain & (bitBoardEmpty[uppermostStone+4] >> 4), 4, winningPosition);
			return true;
		}
		// Diagonally / at the top (basically shift and then check vertically).
		domain = bitBoard[uppermostStone] & (bitBoard[uppermostStone+1] << 1) & (bitBoard[uppermostStone+2] << 2) & (bitBoard[uppermostStone+3] << 3) & (bitBoardEmpty[uppermostStone+4] << 4);
		if(domain != 0) {
			AI::setwinningPosition(leftColumnLimit, rightColumnLimit, uppermostStone+4, domain & (bitBoardEmpty[uppermostStone+4] << 4), -4, winningPosition);
			return true;
		}
		// Vertically at the bottom.
		domain = bitBoardEmpty[lowermostStone-4] & bitBoard[lowermostStone-3] & bitBoard[lowermostStone-2] & bitBoard[lowermostStone-1] & bitBoard[lowermostStone];
		if(domain != 0) {
			AI::setwinningPosition(leftColumnLimit, rightColumnLimit, lowermostStone-4, domain & bitBoardEmpty[lowermostStone-4], 0, winningPosition);
			return true;
		}
		// Diagonally \ at the bottom (basically shift and then check vertically).
		domain = (bitBoardEmpty[lowermostStone-4] << 1) & bitBoard[lowermostStone-3] & (bitBoard[lowermostStone-2] >> 1) & (bitBoard[lowermostStone-1] >> 2) & (bitBoard[lowermostStone] >> 3);
		if(domain != 0) {
			AI::setwinningPosition(leftColumnLimit, rightColumnLimit, lowermostStone-4, domain & (bitBoardEmpty[lowermostStone-4] << 1), -1, winningPosition);
			return true;
		}
		// Diagonally / at the bottom (basically shift and then check vertically).
		domain = (bitBoardEmpty[lowermostStone-4] >> 1) & bitBoard[lowermostStone-3] & (bitBoard[lowermostStone-2] << 1) & (bitBoard[lowermostStone-1] << 2) & (bitBoard[lowermostStone] << 3);
		if(domain != 0) {
			AI::setwinningPosition(leftColumnLimit, rightColumnLimit, lowermostStone-4, domain & (bitBoardEmpty[lowermostStone-4] >> 1), 1, winningPosition);
			return true;
		}
	}
	// Horizontally.
	for(int i = uppermostStone; i <= lowermostStone; i++) {
		domain = ((bitBoardEmpty[i] << 1) | (bitBoardEmpty[i] >> 4)) & bitBoard[i] & (bitBoard[i] >> 1) & (bitBoard[i] >> 2) & (bitBoard[i] >> 3);
		if(domain != 0) {
			if(!AI::setwinningPosition(leftColumnLimit, rightColumnLimit, i, domain & (bitBoardEmpty[i] << 1), -1, winningPosition)) {
				AI::setwinningPosition(leftColumnLimit, rightColumnLimit, i, domain & (bitBoardEmpty[i] >> 4), 4, winningPosition);				
			}
			return true;
		}
	}
	return false;
}

int AI::numberOfFoursOpen(Board *board, std::vector<int> &bitBoard, std::vector<int> &bitBoardEmpty) {
	int numberOfDomains = 0;
	int domain;
	int uppermostStone = board->searchRectangle.upperLeft.y;
	int lowermostStone = board->searchRectangle.lowerRight.y;	
	for(int i = uppermostStone+1; i <= lowermostStone-4; i++) {
		// Vertically.
		domain = (bitBoardEmpty[i-1] & bitBoardEmpty[i+4]) & bitBoard[i] & bitBoard[i+1] & bitBoard[i+2] & bitBoard[i+3];
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally \ (basically shift and then check vertically).
		domain = ((bitBoardEmpty[i-1] << 1) & (bitBoardEmpty[i+4] >> 4)) & bitBoard[i] & (bitBoard[i+1] >> 1) & (bitBoard[i+2] >> 2) & (bitBoard[i+3] >> 3);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally / (basically shift and then check vertically).
		domain = ((bitBoardEmpty[i-1] >> 1) & (bitBoardEmpty[i+4] << 4)) & bitBoard[i] & (bitBoard[i+1] << 1) & (bitBoard[i+2] << 2) & (bitBoard[i+3] << 3);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
	}
	// Horizontally.
	for(int i = uppermostStone; i <= lowermostStone; i++) {
		domain = ((bitBoardEmpty[i] << 1) & (bitBoardEmpty[i] >> 4)) & bitBoard[i] & (bitBoard[i] >> 1) & (bitBoard[i] >> 2) & (bitBoard[i] >> 3);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
	}
	return numberOfDomains;
}

int AI::numberOfFoursHalfOpen(Board *board, std::vector<int> &bitBoard, std::vector<int> &bitBoardEmpty) {
	int numberOfDomains = 0;
	int domain;
	int uppermostStone = board->searchRectangle.upperLeft.y;
	int lowermostStone = board->searchRectangle.lowerRight.y;
	for(int i = uppermostStone+1; i <= lowermostStone-4; i++) {
		// Vertically.
		domain = (bitBoardEmpty[i-1] | bitBoardEmpty[i+4]) & bitBoard[i] & bitBoard[i+1] & bitBoard[i+2] & bitBoard[i+3];
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally \ (basically shift and then check vertically).
		domain = ((bitBoardEmpty[i-1] << 1) | (bitBoardEmpty[i+4] >> 4)) & bitBoard[i] & (bitBoard[i+1] >> 1) & (bitBoard[i+2] >> 2) & (bitBoard[i+3] >> 3);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally / (basically shift and then check vertically).
		domain = ((bitBoardEmpty[i-1] >> 1) | (bitBoardEmpty[i+4] << 4)) & bitBoard[i] & (bitBoard[i+1] << 1) & (bitBoard[i+2] << 2) & (bitBoard[i+3] << 3);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
	}
	// We have to treat the top and bottom separately.
	if(uppermostStone+4 <= lowermostStone) {
		// Vertically at the top.
		domain = bitBoard[uppermostStone] & bitBoard[uppermostStone+1] & bitBoard[uppermostStone+2] & bitBoard[uppermostStone+3] & bitBoardEmpty[uppermostStone+4];
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally \ at the top (basically shift and then check vertically).
		domain = bitBoard[uppermostStone] & (bitBoard[uppermostStone+1] >> 1) & (bitBoard[uppermostStone+2] >> 2) & (bitBoard[uppermostStone+3] >> 3) & (bitBoardEmpty[uppermostStone+4] >> 4);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally / at the top (basically shift and then check vertically).
		domain = bitBoard[uppermostStone] & (bitBoard[uppermostStone+1] << 1) & (bitBoard[uppermostStone+2] << 2) & (bitBoard[uppermostStone+3] << 3) & (bitBoardEmpty[uppermostStone+4] << 4);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Vertically at the bottom.
		domain = bitBoardEmpty[lowermostStone-4] & bitBoard[lowermostStone-3] & bitBoard[lowermostStone-2] & bitBoard[lowermostStone-1] & bitBoard[lowermostStone];
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally \ at the bottom (basically shift and then check vertically).
		domain = (bitBoardEmpty[lowermostStone-4] << 1) & bitBoard[lowermostStone-3] & (bitBoard[lowermostStone-2] >> 1) & (bitBoard[lowermostStone-1] >> 2) & (bitBoard[lowermostStone] >> 3);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally / at the bottom (basically shift and then check vertically).
		domain = (bitBoardEmpty[lowermostStone-4] >> 1) & bitBoard[lowermostStone-3] & (bitBoard[lowermostStone-2] << 1) & (bitBoard[lowermostStone-1] << 2) & (bitBoard[lowermostStone] << 3);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
	}
	// Horizontally.
	for(int i = uppermostStone; i <= lowermostStone; i++) {
		domain = ((bitBoardEmpty[i] << 1) | (bitBoardEmpty[i] >> 4)) & bitBoard[i] & (bitBoard[i] >> 1) & (bitBoard[i] >> 2) & (bitBoard[i] >> 3);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
	}
	return numberOfDomains;
}

int AI::numberOfThreesOpen(Board *board, std::vector<int> &bitBoard, std::vector<int> &bitBoardEmpty) {
	int numberOfDomains = 0;
	int domain;
	int uppermostStone = board->searchRectangle.upperLeft.y;
	int lowermostStone = board->searchRectangle.lowerRight.y;
	for(int i = uppermostStone+1; i <= lowermostStone-3; i++) {
		// Vertically.
		domain = (bitBoardEmpty[i-1] & bitBoardEmpty[i+3]) & bitBoard[i] & bitBoard[i+1] & bitBoard[i+2];
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally \ (basically shift and then check vertically).
		domain = ((bitBoardEmpty[i-1] << 1) & (bitBoardEmpty[i+3] >> 3)) & bitBoard[i] & (bitBoard[i+1] >> 1) & (bitBoard[i+2] >> 2);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally / (basically shift and then check vertically).
		domain = ((bitBoardEmpty[i-1] >> 1) & (bitBoardEmpty[i+3] << 3)) & bitBoard[i] & (bitBoard[i+1] << 1) & (bitBoard[i+2] << 2);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
	}
	// Horizontally.
	for(int i = uppermostStone; i <= lowermostStone; i++) {
		domain = ((bitBoardEmpty[i] << 1) & (bitBoardEmpty[i] >> 3)) & bitBoard[i] & (bitBoard[i] >> 1) & (bitBoard[i] >> 2);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
	}
	return numberOfDomains;
}

int AI::numberOfThreesHalfOpen(Board *board, std::vector<int> &bitBoard, std::vector<int> &bitBoardEmpty) {
	int numberOfDomains = 0;
	int domain;
	int uppermostStone = board->searchRectangle.upperLeft.y;
	int lowermostStone = board->searchRectangle.lowerRight.y;
	for(int i = uppermostStone+1; i <= lowermostStone-3; i++) {
		// Vertically.
		domain = (bitBoardEmpty[i-1] | bitBoardEmpty[i+3]) & bitBoard[i] & bitBoard[i+1] & bitBoard[i+2];
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally \ (basically shift and then check vertically).
		domain = ((bitBoardEmpty[i-1] << 1) | (bitBoardEmpty[i+3] >> 3)) & bitBoard[i] & (bitBoard[i+1] >> 1) & (bitBoard[i+2] >> 2);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally / (basically shift and then check vertically).
		domain = ((bitBoardEmpty[i-1] >> 1) | (bitBoardEmpty[i+3] << 3)) & bitBoard[i] & (bitBoard[i+1] << 1) & (bitBoard[i+2] << 2);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
	}
	// We have to treat the top and bottom separately.
	if(uppermostStone+3 <= lowermostStone) {
		// Vertically at the top.
		domain = bitBoard[uppermostStone] & bitBoard[uppermostStone+1] & bitBoard[uppermostStone+2] & bitBoardEmpty[uppermostStone+3];
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally \ at the top (basically shift and then check vertically).
		domain = bitBoard[uppermostStone] & (bitBoard[uppermostStone+1] >> 1) & (bitBoard[uppermostStone+2] >> 2) & (bitBoardEmpty[uppermostStone+3] >> 3);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally / at the top (basically shift and then check vertically).
		domain = bitBoard[uppermostStone] & (bitBoard[uppermostStone+1] << 1) & (bitBoard[uppermostStone+2] << 2) & (bitBoardEmpty[uppermostStone+3] << 3);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Vertically at the bottom.
		domain = bitBoardEmpty[lowermostStone-3] & bitBoard[lowermostStone-2] & bitBoard[lowermostStone-1] & bitBoard[lowermostStone];
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally \ at the bottom (basically shift and then check vertically).
		domain = (bitBoardEmpty[lowermostStone-3] << 1) & bitBoard[lowermostStone-2] & (bitBoard[lowermostStone-1] >> 1) & (bitBoard[lowermostStone] >> 2);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally / at the bottom (basically shift and then check vertically).
		domain = (bitBoardEmpty[lowermostStone-3] >> 1) & bitBoard[lowermostStone-2] & (bitBoard[lowermostStone-1] << 1) & (bitBoard[lowermostStone] << 2);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
	}
	// Horizontally.
	for(int i = uppermostStone; i <= lowermostStone; i++) {
		domain = ((bitBoardEmpty[i] << 1) | (bitBoardEmpty[i] >> 3)) & bitBoard[i] & (bitBoard[i] >> 1) & (bitBoard[i] >> 2);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
	}
	return numberOfDomains;
}

int AI::numberOfTwosOpen(Board *board, std::vector<int> &bitBoard, std::vector<int> &bitBoardEmpty) {
	int numberOfDomains = 0;
	int domain;
	int uppermostStone = board->searchRectangle.upperLeft.y;
	int lowermostStone = board->searchRectangle.lowerRight.y;
	for(int i = uppermostStone+1; i <= lowermostStone-2; i++) {
		// Vertically.
		domain = (bitBoardEmpty[i-1] & bitBoardEmpty[i+2]) & bitBoard[i] & bitBoard[i+1];
		while(domain > 0) {
			domain &= domain-1;
			numberOfDomains++;
		}
		// Diagonally \ (basically shift and then check vertically).
		domain = ((bitBoardEmpty[i-1] << 1) & (bitBoardEmpty[i+2] >> 2)) & bitBoard[i] & (bitBoard[i+1] >> 1);
		while(domain > 0) {
			domain &= domain-1;
			numberOfDomains++;
		}
		// Diagonally / (basically shift and then check vertically).
		domain = ((bitBoardEmpty[i-1] >> 1) & (bitBoardEmpty[i+2] << 2)) & bitBoard[i] & (bitBoard[i+1] << 1);
		while(domain > 0) {
			domain &= domain-1;
			numberOfDomains++;
		}
	}
	// Horizontally.
	for(int i = uppermostStone; i <= lowermostStone; i++) {
		domain = ((bitBoardEmpty[i] << 1) & (bitBoardEmpty[i] >> 2)) & bitBoard[i] & (bitBoard[i] >> 1);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
	}
	return numberOfDomains;
}

int AI::numberOfTwosHalfOpen(Board *board, std::vector<int> &bitBoard, std::vector<int> &bitBoardEmpty) {
	int numberOfDomains = 0;
	int domain;
	int uppermostStone = board->searchRectangle.upperLeft.y;
	int lowermostStone = board->searchRectangle.lowerRight.y;
	for(int i = uppermostStone+1; i <= lowermostStone-2; i++) {
		// Vertically.
		domain = (bitBoardEmpty[i-1] | bitBoardEmpty[i+2]) & bitBoard[i] & bitBoard[i+1];
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally \ (basically shift and then check vertically).
		domain = ((bitBoardEmpty[i-1] << 1) | (bitBoardEmpty[i+2] >> 2)) & bitBoard[i] & (bitBoard[i+1] >> 1);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally / (basically shift and then check vertically).
		domain = ((bitBoardEmpty[i-1] >> 1) | (bitBoardEmpty[i+2] << 2)) & bitBoard[i] & (bitBoard[i+1] << 1);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
	}
	// We have to treat the top and bottom separately.
	if(uppermostStone+2 <= lowermostStone) {
		// Vertically at the top.
		domain = bitBoard[uppermostStone] & bitBoard[uppermostStone+1] & bitBoardEmpty[uppermostStone+2];
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally \ at the top (basically shift and then check vertically).
		domain = bitBoard[uppermostStone] & (bitBoard[uppermostStone+1] >> 1) & (bitBoardEmpty[uppermostStone+2] >> 2);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally / at the top (basically shift and then check vertically).
		domain = bitBoard[uppermostStone] & (bitBoard[uppermostStone+1] << 1) & (bitBoardEmpty[uppermostStone+2] << 2);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Vertically at the bottom.
		domain = bitBoardEmpty[lowermostStone-2] & bitBoard[lowermostStone-1] & bitBoard[lowermostStone];
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally \ at the bottom (basically shift and then check vertically).
		domain = (bitBoardEmpty[lowermostStone-2] << 1) & bitBoard[lowermostStone-1] & (bitBoard[lowermostStone] >> 1);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
		// Diagonally / at the bottom (basically shift and then check vertically).
		domain = (bitBoardEmpty[lowermostStone-2] >> 1) & bitBoard[lowermostStone-1] & (bitBoard[lowermostStone] << 1);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
	}
	// Horizontally.
	for(int i = uppermostStone; i <= lowermostStone; i++) {
		domain = ((bitBoardEmpty[i] << 1) | (bitBoardEmpty[i] >> 2)) & bitBoard[i] & (bitBoard[i] >> 1);
		while(domain > 0) {
			domain &= domain - 1;
			numberOfDomains++;
		}
	}
	return numberOfDomains;
}
